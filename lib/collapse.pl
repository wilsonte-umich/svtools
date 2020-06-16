use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $tmpH
            $inputDir $tmpDir $outputDir $plotDir
            %scnPosCols %scnCnvCols);
my %p = %scnPosCols;
my %c = %scnCnvCols;
my (@globs, $group, $maxMem, $minCnvDelta);
my ($cnvH, $regH);
my (%cnvPos, %cnvCall, %samples, @samples,
	%posFiles, %posData, @pos, %revPos);
my ($regId, $maxEnd, $regStart, $wrkChrom, $genomeBinSize) = (0, 0);

# manage options
sub setOptions_collapse {
    setOptionValue(\$group,      'group');
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);   
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$minCnvDelta,'min-cnv-delta', 5);
    #-------------------------------------
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
}

# main execution block
sub svtools_collapse {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";
    $globs[1] or die "$error: expected two or more svtools.crunch.cnvs files as input\n";
    
    # initialize
    setOptions_collapse();    
    print STDERR "$utility collapse: " . getTime(), "\n";
	print STDERR "comparing CNVs called by scan/crunch between samples\n";    
	foreach my $glob(@globs){
		foreach my $cnvFile(glob($glob)){
			$cnvFile =~ m/(.*svtools)\.crunch\.cnvs\.(.+)\.bgz/ or die "not an svtools.crunch.cnvs file: $cnvFile\n";
            $samples{$2}++;
			$posFiles{$2} = "$1.scan.positions.$2.bgz";
		}
	}
    @samples = keys %samples;
    print STDERR "  ", join("\n  ", @samples), "\n";

	# open streams
	my $regFile = getOutFile("regions", $group, "bgz");     
	openOutputStream($regFile, \$regH, undef, undef,
				     "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
	
	############
	openInputStream("gunzip -c @globs", \$cnvH, undef, "awk '\$1==\"chr19\"' | sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n");

	# do the work
	while (my $line = <$cnvH>) {
		chomp $line;
		my @f = split("\t", $line);
		if (($wrkChrom and $wrkChrom ne $f[$c{chrom}]) or 
			($maxEnd   and $maxEnd    < $f[$c{start}]) ) {		
			processRegion();
		}		
		if (!$wrkChrom or $wrkChrom ne $f[$c{chrom}]) {
			$wrkChrom = $f[$c{chrom}];
			loadChromPos($wrkChrom);
		}		
		defined $regStart or $regStart = $f[$c{start}];
		$maxEnd >= $f[$c{end}] or $maxEnd = $f[$c{end}];
        my $smp = $f[$c{sample}];
		push @{$cnvPos{$smp}}, \@f;
        abs($f[$c{deltaCDF_sgn}]) >= $minCnvDelta and push @{$cnvCall{$smp}}, \@f;
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
            my ($maxCnvDltIdx,   $maxPosDltIdx,
				$maxCnvDltOther, $maxPosDltOther,
				$isGap, @call) = (0, 0, 0, 0, 0, 0);
            foreach my $cnv(@{$cnvPos{$idxSmp}}){ # get the best CNV span for this sample
                my $dlt = abs($$cnv[$c{deltaCDF_sgn}]);
                $maxCnvDltIdx >= $dlt or $maxCnvDltIdx = $dlt;
                if ($maxCnvDltIdx == $dlt) { # take most extreme span as the cnv call for this region for this sample
                    @call = map { $$cnv[$c{$_}] } qw(
                        end nBins nFrags
                        fracLowQual tLens isGap n_log_p deltaCDF_abs deltaCDF_sgn
                        deltaCDF_alt cnvType cnvSize cnvFrac
                    );
                }                
            }			
            foreach my $othSmp(@samples){ # get the best CNV span for all other samples
                $othSmp eq $idxSmp and next;
                $cnvPos{$othSmp} or next;
                foreach my $cnv(@{$cnvPos{$othSmp}}){
                    my $dlt = abs($$cnv[$c{deltaCDF_sgn}]);
                    $maxCnvDltOther >= $dlt or $maxCnvDltOther = $dlt;  
                }                
            }			
			my %gapSamples;
			for (my $pos=$regStart+1; $pos<=$maxEnd; $pos+=$genomeBinSize){ # get the max position delta and gap information for this sample
				my $i = $revPos{$pos};
				defined $i or die "missing index for position: $wrkChrom, $pos\n";
				my $idxData = $posData{$idxSmp}[$i];
				my $dlt = ($idxData ? abs($$idxData[1]) : 0);
				$maxPosDltIdx >= $dlt or $maxPosDltIdx = $dlt;
				$isGap = ($isGap || ($idxData ? $$idxData[0] : 1));					
				foreach my $othSmp(@samples){ # get the best CNV span for all other samples
					$othSmp eq $idxSmp and next;
					my $othData = $posData{$othSmp}[$i];
					
					###############
					#print $othData, "\n";
					
					my $dlt = ($othData ? abs($$othData[1]) : 0);				
					$maxPosDltOther >= $dlt or $maxPosDltOther = $dlt;
					$isGap = ($othData ? $$othData[0] : 1);
					$isGap and $gapSamples{$othSmp}++;					
				}				
			}
			my $nOtherGap = scalar(keys %gapSamples);
            print $regH join("\t", $cnvCall{$idxSmp}[0][$c{chrom}], $regStart, $maxEnd, $regId, 1+($maxEnd-$regStart-1)/$genomeBinSize, '.',
                                   $nPos, $nCall,
                                   $idxSmp, $maxCnvDltIdx, $maxPosDltIdx, $maxCnvDltOther, $maxPosDltOther,
								   $isGap, $nOtherGap,
                                   @call), "\n";		            
        }  
    }
	$regStart = undef;	
	$maxEnd = 0;
    %cnvPos = %cnvCall = ();
}

sub loadChromPos {
	my ($chrom) = @_;
	print STDERR "  $chrom loading ";
	%posData = @pos = ();
	foreach my $sample(@samples){
		print STDERR ".";
		openInputStream("tabix $posFiles{$sample} $chrom", \my $tmpH);
		my $i = 0;
		while (my $line = <$tmpH>) {
			chomp $line;
			my @f = split("\t", $line);
			if (defined $pos[$i] and $pos[$i] != $f[$p{end}]) {
				die "mismatched 'end' column at chrom=$chrom, i=$i, sample=$sample:\n".
				    "found $f[$p{end}], expected $pos[$i]\n";
			}
			defined $pos[$i] or $pos[$i] = $f[$p{end}];			
			$posData{$sample}[$i] = [map {$f[$p{$_}]} qw(isGap deltaCDF_sgn)];
			$i++;
		}	
	}
	%revPos = map { $pos[$_] => $_ } 0..$#pos;
    $genomeBinSize = $pos[1] - $pos[0];	
	print STDERR " working\n";
}



## determine the genome bin size from the first file (must be same for all files)
#sub setGenomeBinSize {
#	openInputStream("gunzip -c @globs", \my $inH);
#	my $line1 = <$inH>;
#	my $line2 = <$inH>;	
#	closeHandles($inH);	
#	my @line1 = split("\t", $line1);
#	my @line2 = split("\t", $line2);
#	$genomeBinSize = $line2[$c0{end}] - $line1[$c0{end}];
#	print STDERR "genomeBinSize = $genomeBinSize\n";
#}

## assemble the command to collapse bins into arrays of non-gap samples
#sub setCollapseStream {
#	my ($chrom) = @_;
#	my $groupCols = join(",", map { $c1{$_} } qw(chrom end));
#	@cllpsCols = qw(sample nFrags deltaCDF_ref deltaCDF_alt cnvSize cnvFrac cnvType cnvType2);
#    my $i = 2;
#	%cc = map { $_ => $i++ } @cllpsCols;
#	my $cllpsCols = join(",", map { $c1{$_} } @cllpsCols);
#	my $cllpsType = join(",", ('collapse') x @cllpsCols);
#    $clpStream =
#		"gunzip -c @globs | ".
#		"awk '\$1==\"$chrom\"' | ".
#		"grep -v gap | ". # ignore bins with insufficient data for a sample
#		"sort -T $tmpDir -S $maxMem"."b -k$c1{end},$c1{end}n | ". 
#		"groupBy -g $groupCols -c $cllpsCols -o $cllpsType | ". # execute the collapse
#		"grep 'gain\\|loss'";	# ignore bins that were 'ref' for all samples (could still have 'del,ref' and similar)
#}    

## thread the collapse bins data to find contiguous runs of CNV bins (in any sample)
#sub analyzeBins {
#	$prevF = undef;
#	while (my $line = <$clpH>) {
#		chomp $line;
#		my @f = split("\t", $line);
#		if ($prevF and # find breakpoints between CNV spans
#			$$prevF[1] != $f[1] - $genomeBinSize) {		
#			processBinGroup();
#		}
#		defined $cnvStart or $cnvStart = $f[1];
#		$prevF = \@f;		
#		my %d = map { $_ => [split(",", $f[$cc{$_}])] } keys %cc; # break the collapse into arrays
#		my $i = 0;
#		$d{smpIs} = { map { $_ => $i++ } @{$d{sample}} }; # set the array indices per sample
#		push @bins, \%d; # save the split bin data
#		foreach my $smp(@{$d{sample}}){ # record how many samples had CNV calls in the current span (not bin)
#			my $smpI = $d{smpIs}{$smp};
#			$d{cnvType}[$smpI] ne 'ref' and $cnvSmp{$smp}++; # and count how many bins were positive for each sample
#		}
#	}
#	processBinGroup();
#}

## print a line for every called sample for a contiguous run of CNV bins
#sub processBinGroup {
#	my $nBins = @bins;
#	if ($nBins >= $minCnvBins) {
#		my ($nPosSmp, $nCllSmp, %toPrint) = (scalar(keys %cnvSmp), 0);
#		foreach my $idxSmp(keys %cnvSmp){
#			$cnvSmp{$idxSmp} >= $minCnvBins or next;
#			$nCllSmp++; # sample had sufficient cnv bins, call it, noting that bins may not be contiguous		
#            my ($maxDltIdx, $maxDltOther, @call) = (-1000, -1000);            
#			foreach my $bin(@bins){
#				my $idxSmpI = $$bin{smpIs}{$idxSmp};
#				if (defined $idxSmpI) { # find the bins for which this sample had data
#					my $dlt = $$bin{deltaCDF_ref}[$idxSmpI];
#					$maxDltIdx >= $dlt or $maxDltIdx = $dlt; # find the most extreme deltaMle for this sample within region
#					if ($maxDltIdx == $dlt) { # take bin with most extreme value as the cnv call for this region for this sample
#						@call = map { $$bin{$_}[$idxSmpI] } qw(nFrags deltaCDF_alt cnvSize cnvFrac cnvType cnvType2);
#					}
#				}            
#				foreach my $smp(@{$$bin{sample}}){
#					$smp eq $idxSmp and next;      # find the most extreme KS p-value for all other samples within region
#					my $smpI = $$bin{smpIs}{$smp}; # whether or not they were called as CNVs
#					my $dlt = $$bin{deltaCDF_ref}[$smpI];
#					$maxDltOther >= $dlt or $maxDltOther = $dlt;
#				}	
#			}
#			$toPrint{$idxSmp} = join("\t", $idxSmp, $cnvSmp{$idxSmp}, $maxDltIdx, $maxDltOther, @call);
#		}
#		my $fracGapSmp = 0; # determine the average fraction of gap samples over all bins in region
#		foreach my $bin(@bins){ $fracGapSmp += ($totSmp - @{$$bin{sample}}) / $totSmp }
#		$fracGapSmp /= $nBins;
#		foreach my $idxSmp(keys %toPrint){ # commit the results for every called CNV in region
#			$cnvId++;
#			print $regH join("\t", $$prevF[0], $cnvStart-1, $$prevF[1], $cnvId, $nBins, '.',
#					         $nPosSmp, $nCllSmp, $fracGapSmp,
#							 $toPrint{$idxSmp}), "\n";	
#		}
#	}
#	$cnvStart = undef; # reset for next region
#	@bins 	  = ();
#	%cnvSmp   = ();
#}

1;

#sub processBinGroup {
#	my $nBins = @bins;
#	if ($nBins >= $minCnvBins) {
#		my ($nPosSmp, $nCllSmp, %toPrint) = (scalar(keys %cnvSmp), 0);
#		foreach my $idxSmp(keys %cnvSmp){
#			$cnvSmp{$idxSmp} >= $minCnvBins or next;
#			$nCllSmp++; # sample had sufficient cnv bins, call it, noting that bins may not be contiguous
#            #my ($maxPIdx, $maxPOther, @call) = (0, 0);			
#            my ($maxDltIdx, $maxDltOther, @call) = (-1000, -1000);            
#
#			foreach my $bin(@bins){
#				my $idxSmpI = $$bin{smpIs}{$idxSmp};
#				if (defined $idxSmpI) { # find the bins for which this sample had data
#            
#
#					my $p = $$bin{nlog_ks_p_ref}[$idxSmpI];
#					$maxPIdx >= $p or $maxPIdx = $p; # find the most extreme KS p-value for this sample within region
#					if ($maxPIdx == $p) { # take bin with most extreme p as the cnv call for this region for this sample
#						@call = map { $$bin{$_}[$idxSmpI] } qw(nlog_ks_p_alt nFrags dupSize nDupFrag delInsSize nDelInsFrag);
#					}
#				}
#				foreach my $smp(@{$$bin{sample}}){
#					$smp eq $idxSmp and next;      # find the most extreme KS p-value for all other samples within region
#					my $smpI = $$bin{smpIs}{$smp}; # whether or not they were called as CNVs
#					my $p = $$bin{nlog_ks_p_ref}[$smpI];
#					$maxPOther >= $p or $maxPOther = $p;
#				}	
#			}
#			$toPrint{$idxSmp} = join("\t", $idxSmp, $cnvSmp{$idxSmp}, $maxPIdx, $maxPOther, @call);
#		}
#		my $fracGapSmp = 0; # determine the average fraction of gap samples over all bins in region
#		foreach my $bin(@bins){ $fracGapSmp += ($totSmp - @{$$bin{sample}}) / $totSmp }
#		$fracGapSmp /= $nBins;
#		foreach my $idxSmp(keys %toPrint){ # commit the results for every called CNV in region
#			$cnvId++;
#			print $regH join("\t", $$prevF[0], $cnvStart-1, $$prevF[1], $cnvId, $nBins, '.',
#					         $nPosSmp, $nCllSmp, $fracGapSmp,
#							 $toPrint{$idxSmp}), "\n";	
#		}
#	}
#	$cnvStart = undef; # reset for next region
#	@bins 	  = ();
#	%cnvSmp   = ();
#}
