use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $tmpH
            $inputDir $tmpDir $outputDir $plotDir
            %biasSnpCols %biasCnvCols);
my %p = %biasSnpCols;
my %c = %biasCnvCols;
my (@globs, $group, $maxMem, $minPRatio);
my ($cnvH, $regH);
my (%cnvPos, %cnvCall, %samples, @samples,
	%snpFiles, %snpData, @snps, %revSnps);
my ($regId, $maxEnd, $regStart, $wrkChrom) = (0, 0);

# manage options
sub setOptions_collapse_bias {
    setOptionValue(\$group,      'group');
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);   
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$minPRatio,  'min-pRatio',   10);
    #-------------------------------------
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
}

# main execution block
sub svtools_collapse_bias {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";
    $globs[1] or die "$error: expected two or more svtools.bias.cnvs files as input\n";
    
    # initialize
    setOptions_collapse_bias();    
    print STDERR "$utility collapse_bias: " . getTime(), "\n";
	print STDERR "comparing CNVs called by genotype/bias between samples\n";    
	foreach my $glob(@globs){
		foreach my $cnvFile(glob($glob)){
			$cnvFile =~ m/(.*svtools)\.bias\.cnvs\.(.+)\.bgz/ or die "not an svtools.bias.cnvs file: $cnvFile\n";
            $samples{$2}++;
			$snpFiles{$2} = "$1.bias.snps.$2.bgz";
		}
	}
    @samples = keys %samples;
    print STDERR "  ", join("\n  ", @samples), "\n";

	# open streams
	my $regFile = getOutFile("regions", $group, "bgz");     
	openOutputStream($regFile, \$regH, undef, undef,
				     "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
	openInputStream("gunzip -c @globs", \$cnvH, undef, "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n");

	# do the work
	while (my $line = <$cnvH>) {
		chomp $line;
		my @f = split("\t", $line);
		if (($wrkChrom and $wrkChrom ne $f[$c{chrom}]) or 
			($maxEnd   and $maxEnd    < $f[$c{end}]) ) {		
			processRegion();
		}		
		if (!$wrkChrom or $wrkChrom ne $f[$c{chrom}]) {
			$wrkChrom = $f[$c{chrom}];
			loadChromSnps($wrkChrom);
		}		
		defined $regStart or $regStart = $f[$c{start}];
		$maxEnd >= $f[$c{end}] or $maxEnd = $f[$c{end}];
        my $smp = $f[$c{sample}];
		push @{$cnvPos{$smp}}, \@f;
        $f[$c{pRatio}] >= $minPRatio and push @{$cnvCall{$smp}}, \@f;
	}
	processRegion();

    # finish up
	closeHandles($cnvH, $regH);
	system("tabix -p bed $regFile");
}

sub processRegion {
	my $nCall = scalar(keys %cnvCall);
	if($nCall){
        my $nPos = scalar(keys %cnvPos);
        $regId++;
        foreach my $idxSmp(keys %cnvCall){        
            my ($maxPRatioIdx, $maxPRatioOther, $cnvSize, @call) = (0, 0, 0);
			my ($sign, $minCnvI, $maxCnvI, $maxOtherBias) = (0, 0, 0, 0);
            foreach my $cnv(@{$cnvPos{$idxSmp}}){ # get the best CNV span for this sample
                my $pRat = $$cnv[$c{pRatio}];
                $maxPRatioIdx >= $pRat or $maxPRatioIdx = $pRat;
                if ($maxPRatioIdx == $pRat) { # take most extreme span as the cnv call for this region for this sample
                    @call = map { $$cnv[$c{$_}] } qw(
                        start end cnvId nSnps nRef nAlt biasType  
                    );
					my $bias = ($$cnv[$c{nRef}] - $$cnv[$c{nAlt}]) /
							   ($$cnv[$c{nRef}] + $$cnv[$c{nAlt}]);
					$sign = $bias < 0 ? -1 : 1;
					push @call, ($$cnv[$c{end}] - $$cnv[$c{start}],
								 int(abs($bias) / 0.001 + 0.5) * 0.001 * $sign);
					$minCnvI = $revSnps{$$cnv[$c{start}] + 1};
					$maxCnvI = $revSnps{$$cnv[$c{end}]};
                }                
            }
            foreach my $othSmp(@samples){ # get the best CNV span for all other samples
                $othSmp eq $idxSmp and next;
                if($cnvPos{$othSmp}){
					foreach my $cnv(@{$cnvPos{$othSmp}}){
						my $pRat = $$cnv[$c{pRatio}];
						$maxPRatioOther >= $pRat or $maxPRatioOther = $pRat;  
					}
				}
				my ($nRef, $nAlt) = (0, 0);
				foreach my $i($minCnvI..$maxCnvI){
					my $snp = $snpData{$othSmp}[$i];
					$nRef += $$snp[0];
					$nAlt += $$snp[1];
				}
				my $nReads = $nRef + $nAlt;
				my $otherBias = $nReads ? (($nRef - $nAlt) / $nReads) : 0;
				my $otherSign = $otherBias < 0 ? -1 : 1;
				if ($sign == $otherSign) {
					$otherBias = int(abs($otherBias) / 0.001 + 0.5) * 0.001;
				} else {
					$otherBias = 0;
				}
				$maxOtherBias >= $otherBias or $maxOtherBias = $otherBias;
            }
            print $regH join("\t", $cnvCall{$idxSmp}[0][$c{chrom}], $regStart, $maxEnd, $regId, $maxCnvI - $minCnvI + 1, '.',
                                   $nPos, $nCall,
                                   $idxSmp, int($maxPRatioIdx/0.001+0.5)*0.001, int($maxPRatioOther/0.001+0.5)*0.001,
                                   @call, $maxOtherBias * $sign), "\n";  
        }  
    }
	$regStart = undef;	
	$maxEnd = 0;
    %cnvPos = %cnvCall = ();
}

sub loadChromSnps {
	my ($chrom) = @_;
	print STDERR "  $chrom\tloading ";
	%snpData = @snps = ();
	foreach my $sample(@samples){
		print STDERR ".";
		openInputStream("tabix $snpFiles{$sample} $chrom", \my $tmpH);
		my $i = 0;
		while (my $line = <$tmpH>) {
			chomp $line;
			my @f = split("\t", $line);
			if (defined $snps[$i] and $snps[$i] != $f[$p{end}]) {
				die "mismatched 'end' column at chrom=$chrom, i=$i, sample=$sample:\n".
				    "found $f[$p{end}], expected $snps[$i]\n";
			}
			defined $snps[$i] or $snps[$i] = $f[$p{end}];			
			$snpData{$sample}[$i] = [map {$f[$p{$_}]} qw(nRef nAlt)];
			$i++;
		}	
	}
	%revSnps = map { $snps[$_] => $_ } 0..$#snps;
	print STDERR " working\n";
}
