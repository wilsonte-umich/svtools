use strict;
use warnings;
use List::Util qw(shuffle);

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $tmpH
            $inputDir $tmpDir $outputDir $plotDir);
require "$libDir/exclude_subs.pl";
my (@globs, $sample, $pdfDir, $excludeFile, $maxMem,
    $minQual, $genBinSize, $maxCnvSize,
    $minScanPairs, $maxScanPairs, $minNegPairs,
	$chroms, $region);
my (%chroms);
my ($minTLen, $maxTLen);
my ($tLenH, $scnH, $binH);
my ($libraryN, %libNs) = (0);
my ($bin0, $leftBinI, $rghtBinI, @nTot, @nLowQual, @binTLen, @binLibs_XXX, @enc) = (1);
my ($prevChrom);

# manage options
sub setOptions_scan {
    setOptionValue(\$sample,     'sample');
    setOptionValue(\$pdfDir,     'pdf-dir');  
    setOptionValue(\$excludeFile,'exclude-file');    
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);      
    setOptionValue(\$outputDir,  'output-dir');  
    setOptionValue(\$minQual,    'min-qual',        5);
    setOptionValue(\$genBinSize, 'genome-bin-size', 1000);
    setOptionValue(\$maxCnvSize, 'max-cnv-size',    10000);
    setOptionValue(\$minScanPairs,'min-scan-pairs', 10);
	setOptionValue(\$maxScanPairs,'max-scan-pairs', 100);
	setOptionValue(\$minNegPairs, 'min-neg-pairs',   100);
	setOptionValue(\$chroms,     'chroms');
    setOptionValue(\$region,     'region');
    #-------------------------------------
    -d $pdfDir or die "$error: $pdfDir does not exist or is not a directory\n";
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    $minQual =~ m|\D| and die "$error: min-qual must be an integer MAPQ\n";    
    $genBinSize =~ m|\D| and die "$error: genome-bin-size must be an integer number of bp\n";
    $maxCnvSize =~ m|\D| and die "$error: max-cnv-size must be an integer number of bp\n";
    $minScanPairs =~ m|\D| and die "$error: min-scan-pairs must be an integer number\n";
	$maxScanPairs =~ m|\D| and die "$error: max-scan-pairs must be an integer number\n";
	$minNegPairs =~ m|\D| and die "$error: min-neg-pairs must be an integer number\n";
    #-------------------------------------
    $minTLen = -$maxCnvSize;
    $maxTLen = 2 * $maxCnvSize;
    $maxScanPairs = $maxScanPairs + 1; # e.g. limit 100 leave all above at nFrags = 101
	#-------------------------------------
	($region or $chroms) or die "$error: either region or chroms must be provided\n";
	if ($region) {
		$region =~ m/(chr.+):\d+-\d+/ or die "region must be in form chr:start-end";
		$chroms{$1}++;
	} else {
		%chroms = map { $_ => 1} split(",", $chroms);
	}	
}

# main execution block
sub svtools_scan {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";

    # initialize
    setOptions_scan();    
    print STDERR "$utility scan: " . getTime(), "\n";
    initializeExclude($excludeFile);

    # collect tLens per bin crossing point
	print STDERR "collecting TLENs for each genome bin position\n";
    my $R = $region ? "-R $region" : "";
    my $mrgStream = "samtools merge -u $R - @globs | samtools view -";
    openInputStream($mrgStream, \$scnH);	
    my $tLenFile = getTmpFile("bin_tLens", $sample, "gz");
    openOutputStream($tLenFile, \$tLenH, $TRUE);
    extract_scan();
    closeHandles($scnH, $tLenH);

    # calculating bin deviations
	print STDERR "calculating bin deviations\n";
	my $binFile = getOutFile("positions", $sample, "bgz");
	unless($region){
        openOutputStream($binFile, \$binH, undef, undef,
				"sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    }
	run_scan($tLenFile);
	closeHandles($binH);   
    
    # finish up
    unlink $tLenFile;
	unless($region){
        system("tabix -p bed $binFile");
    }
}

# collect the tLens of all fragments crossing bin points
sub extract_scan {
    my $binBufferSize = 10000;
    while (my $line = <$scnH>) {
		my @aln = split("\t", $line, 12);	
		$aln[11] =~ m/RG:Z:(\S+)/ or die "missing read group in line:\n$line";
        unless($libNs{$1}) {
            $libraryN++;
            $libNs{$1} = $libraryN;
            print STDERR "  library $libraryN = $1\n";
        }        
		my $libN = $libNs{$1} or die "unknown read group '$1' encountered in line:\n$line";
		if(($aln[1] & 0x930) == 0x020 and # forward read of FR pair, not secondary or supplemental
		    $chroms{$aln[2]} and          # on an allowed chromosome
			$aln[6] eq '=' and            # reads on same chromosome		
			$aln[8] >= $minTLen and $aln[8] <= $maxTLen){  # within interrogation size range
                $prevChrom and $aln[2] ne $prevChrom and finishChrom($aln[2]);
                $prevChrom or $prevChrom = $aln[2];	
                my $fPos = $aln[3];
                my $rPos = $fPos + $aln[8];
                my ($leftPos, $rghtPos) = ($fPos <= $rPos) ? ($fPos, $rPos) : ($rPos, $fPos);
                checkExclude($aln[2], $leftPos-1, $rghtPos) and next; # not excluded
                $leftBinI = int($leftPos / $genBinSize) + 1;
                $rghtBinI = int($rghtPos / $genBinSize);
                $leftBinI > $bin0 + $binBufferSize and commitBins($leftBinI - 1);
                my $pairKey = "$libN,$fPos,$aln[8]";
                for (my $binI=$leftBinI; $binI<=$rghtBinI; $binI++){
                    my $i = $binI - $bin0;
					$i < 0 and next;
                    $enc[$i] or $enc[$i] = {};					
                    unless($enc[$i]{$pairKey}){
						$nTot[$i]++;
						if ($aln[4] >= $minQual) { # of sufficient uniqueness
							push @{$binTLen[$i]}, $aln[8];
							#push @{$binLibs[$i]}, $libN;							
						} else {
							$nLowQual[$i]++;
						}
                    }
                    $enc[$i]{$pairKey}++; # suppress duplicate pairs
                }
        }
    }
	finishChrom();
}
sub commitBins {
    my ($maxBin) = @_;
	foreach my $binI($bin0..$maxBin){
		my $i = $binI - $bin0;
		$i < 0 and next;
		my $nFrags = (defined $binTLen[$i] ? @{$binTLen[$i]} : 0);
		my $fracLowQual = $nTot[$i] ? (int((($nLowQual[$i]||0)/$nTot[$i])/0.01+0.5)*0.01) : 0;		
		my ($tLens, $libs);
		if ($nFrags == 0) { # bin has no data (R will ignore)
			($tLens, $libs) = ('NA', 'NA');
		} elsif ($nFrags <= $maxScanPairs) { # typical bin, use as is
			$tLens = join(",", @{$binTLen[$i]});
			#$libs  = join(",", @{$binLibs[$i]});
		} else { # downsample bins with too many pairs
			my @is = (shuffle(0..$nFrags-1))[1..$maxScanPairs];
			$tLens = join(",", @{$binTLen[$i]}[@is]);
			#$libs  = join(",", @{$binLibs[$i]}[@is]);
			$nFrags = $maxScanPairs;
		}
		#print $tLenH join("\t", $prevChrom, $binI, $nFrags, $fracLowQual, $tLens, $libs), "\n";
		print $tLenH join("\t", $prevChrom, $binI, $nFrags, $fracLowQual, $tLens, 'NA'), "\n";
		###################
		#$binI == 143746 and print join("\t", $prevChrom, $binI, $nFrags, $fracLowQual, $tLens, 'NA'), "\n";
	}
	splice(@nTot,    0, $maxBin - $bin0 + 1);	
	splice(@nLowQual,0, $maxBin - $bin0 + 1);	
    splice(@binTLen, 0, $maxBin - $bin0 + 1);
	#splice(@binLibs, 0, $maxBin - $bin0 + 1);
    splice(@enc,     0, $maxBin - $bin0 + 1);
    $bin0 = $maxBin + 1;
}
sub finishChrom {
	my ($chrom) = @_;
	print STDERR "  finishing $prevChrom\n";	
	commitBins($rghtBinI);
	#($bin0, $leftBinI, $rghtBinI, @nTot, @nLowQual, @binTLen, @binLibs, @enc) = (1);
    ($bin0, $leftBinI, $rghtBinI, @nTot, @nLowQual, @binTLen, @enc) = (1);
	$prevChrom = $chrom;
}

# run deviation estimates via R and commit the results
sub run_scan {
	my ($tLenFile) = @_;
    $ENV{sample} = $sample;
	#$ENV{libraries}     = join(",", sort {$libNs{$a} <=> $libNs{$b}} keys %libNs);
    #$ENV{libraryNs}     = join(",", sort {$a <=> $b} values %libNs);    
    $ENV{pdfDir}         = $pdfDir;
    $ENV{tLenFile}       = $tLenFile;
	$ENV{maxCnvSize}     = $maxCnvSize;
    $ENV{min_scan_pairs} = $minScanPairs;
	$ENV{min_neg_pairs}  = $minNegPairs;
    open(my $rH, "-|", "Rscript $libDir/scan.R") or die "could not open scan.R stream: $!\n";
	my ($nStatBins, %n_log_p, %deltaCDF_abs, %deltaCDF_sgn);
    while (my $line = <$rH>) {
        chomp $line;
		my @f = split("\t", $line);
        #my ($chrom, $binI, $nFrags, $fracLowQual, $tLens, 
        #    $isGap, $n_log_p, $deltaCDF_abs, $deltaCDF_sgn)
        my $bin = join("\t",
			$f[0], $f[1]*$genBinSize-1, $f[1]*$genBinSize,
			$sample, $f[2], '.',
			@f[3..8]
		)."\n";
        if ($region) { $f[2] and print STDOUT $bin } else { print $binH $bin }
		if (!$f[5] and $f[3]<0.1) {
			$n_log_p{int($f[6]/0.5+0.5)*0.5}++;
			$deltaCDF_abs{int($f[7]/0.5+0.5)*0.5}++;
			$deltaCDF_sgn{int($f[8]/0.5+0.5)*0.5}++;
			$nStatBins++;
		}
    }
    closeHandles($rH);
	$region and return;
	printScanHist('n_log_p',      \%n_log_p,      $nStatBins);
	printScanHist('deltaCDF_abs', \%deltaCDF_abs, $nStatBins);
	printScanHist('deltaCDF_sgn', \%deltaCDF_sgn, $nStatBins);
}
sub printScanHist {
	my ($name, $hash, $nStatBins) = @_;
	print STDERR "\n$name\tfreq\n";
	foreach my $statVal(sort {$a <=> $b} keys %$hash){
		my $freq = $$hash{$statVal} / $nStatBins;
		print STDERR "$statVal\t$freq\n";
	}
}
            
#Col	Field	Description
#1	QNAME	Query template/pair NAME
#2	FLAG	bitwise FLAG
#3	RNAME	Reference sequence NAME
#4	POS	1-based leftmost POSition/coordinate of clipped sequence
#5	MAPQ	MAPping Quality (Phred-scaled)
#6	CIAGR	extended CIGAR string
#7	MRNM	Mate Reference sequence NaMe (?=? if same as RNAME)
#8	MPOS	1-based Mate POSistion
#9	TLEN	inferred Template LENgth (insert size)
#10	SEQ	query SEQuence on the same strand as the reference
#11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
#12+	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE
#Each bit in the FLAG field is defined as:
#
#Flag	Chr	Description
#0x001	p	the read is paired in sequencing
#0x002	P	the read is mapped in a proper pair
#0x004	u	the query sequence itself is unmapped
#0x008	U	the mate is unmapped
#0x010	r	strand of the query (1 for reverse)
#0x020	R	strand of the mate
#0x040	1	the read is the first read in a pair
#0x080	2	the read is the second read in a pair
#0x100	s	the alignment is not primary
#0x200	f	the read fails platform/vendor quality checks
#0x400	d	the read is either a PCR or an optical duplicate
#0x800      supplementary alignment

1;
