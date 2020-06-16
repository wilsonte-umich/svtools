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
    $minQual, $genBinSize, $cnvSizeStep, $minCnvSize, $maxCnvSize,
    $minMlePairs, $maxMlePairs, $minCnvBins, $minCnvPairs, $region);
my ($minTLen, $maxTLen);
my ($tLenH, $mleH, $binH, $cnvH);
my ($libraryN, %libNs) = (0);
my ($bin0, $leftBinI, $rghtBinI, @binTLen, @binLibs_XXX, @enc) = (1);
my ($prevChrom, $prevBin, @cnvBins);

# manage options
sub setOptions_mle {
    setOptionValue(\$sample,     'sample');
    setOptionValue(\$pdfDir,     'pdf-dir');  
    setOptionValue(\$excludeFile,'exclude-file');    
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);      
    setOptionValue(\$outputDir,  'output-dir');  
    setOptionValue(\$minQual,    'min-qual',     5);
    setOptionValue(\$genBinSize, 'genome-bin-size', 1000);
    setOptionValue(\$cnvSizeStep,'cnv-size-step', 500);
    setOptionValue(\$minCnvSize, 'min-cnv-size', 2000);
    setOptionValue(\$maxCnvSize, 'max-cnv-size', 10000);
    setOptionValue(\$minMlePairs,'min-mle-pairs', 10);
	setOptionValue(\$maxMlePairs,'max-mle-pairs', 100);
	setOptionValue(\$minCnvBins, 'min-cnv-bins', 3);
	setOptionValue(\$minCnvPairs,'min-cnv-pairs', 4);
    setOptionValue(\$region, 'region');
    #-------------------------------------
    -d $pdfDir or die "$error: $pdfDir does not exist or is not a directory\n";
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    $minQual =~ m|\D| and die "$error: min-qual must be an integer MAPQ\n";    
    $genBinSize =~ m|\D| and die "$error: genome-bin-size must be an integer number of bp\n";
    $cnvSizeStep =~ m|\D| and die "$error: cnv-size-step must be an integer number of bp\n";
    $minCnvSize =~ m|\D| and die "$error: min-cnv-size must be an integer number of bp\n";
    $maxCnvSize =~ m|\D| and die "$error: max-cnv-size must be an integer number of bp\n";
    $minCnvSize % $cnvSizeStep and die "$error: min-cnv-size must be a multiple of cnv-size-step\n";
    $maxCnvSize % $cnvSizeStep and die "$error: max-cnv-size must be a multiple of cnv-size-step\n";
    $minMlePairs =~ m|\D| and die "$error: min-mle-pairs must be an integer number\n";
	$maxMlePairs =~ m|\D| and die "$error: max-mle-pairs must be an integer number\n";
	$minCnvBins =~ m|\D| and die "$error: min-cnv-bins must be an integer number\n";
	$minCnvPairs =~ m|\D| and die "$error: min-cnv-pairs must be an integer number\n";
    $region and ($region =~ m/chr.+:\d+-\d+/ or die "region must be in form chr:start-end");
    #-------------------------------------
    $minTLen = -$maxCnvSize;
    $maxTLen = 2 * $maxCnvSize;
    $maxMlePairs = $maxMlePairs + 1; # e.g. limit 100 leave all above at nFrags = 101
}

# main execution block
sub svtools_mle {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";

    # initialize
    setOptions_mle();    
    print STDERR "$utility mle: " . getTime(), "\n";
    initializeExclude($excludeFile);

    # collect tLens per bin crossing point
	print STDERR "collecting TLENs for each genome bin position\n";
    my $R = $region ? "-R $region" : "";
    
    ################
    #$R = "-R chr1";
    
    my $mrgStream = "samtools merge -u $R - @globs | samtools view -";
    openInputStream($mrgStream, \$mleH);	
    my $tLenFile = getTmpFile("bin_tLens", $sample, "gz");
    openOutputStream($tLenFile, \$tLenH, $TRUE);
    extract_mle();
    closeHandles($mleH, $tLenH);

    # perform maximum likelihood analysis
	print STDERR "running MLE analysis\n";
	my $binFile = getOutFile("positions", $sample, "bgz");
	my $cnvFile = getOutFile("cnvs", $sample, "bgz"); 
	unless($region){
        openOutputStream($binFile, \$binH, undef, undef,
				"sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
        openOutputStream($cnvFile, \$cnvH, undef, undef,
                "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    }
	runMle($tLenFile);
	closeHandles($binH, $cnvH);   
    
    # finish up
    unlink $tLenFile;
	unless($region){
        system("tabix -p bed $binFile");
        system("tabix -p bed $cnvFile");
    }
}

# collect the tLens of all fragments crossing bin points
sub extract_mle {
    my $binBufferSize = 10000;
    while (my $line = <$mleH>) {
		my @aln = split("\t", $line, 12);	
		$aln[11] =~ m/RG:Z:(\S+)/ or die "missing read group in line:\n$line";
        unless($libNs{$1}) {
            $libraryN++;
            $libNs{$1} = $libraryN;
            print STDERR "  library $libraryN = $1\n";
        }        
		my $libN = $libNs{$1} or die "unknown read group '$1' encountered in line:\n$line";
		if(($aln[1] & 0x930) == 0x020 and # forward read of FR pair, not secondary or supplemental
			$aln[6] eq '=' and            # on same chromosome		
			$aln[8] >= $minTLen and $aln[8] <= $maxTLen and # within interrogation size range
			$aln[4] >= $minQual){         # of sufficient uniqueness
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
                    $enc[$i] or $enc[$i] = {};
                    unless($enc[$i]{$pairKey}){
                        push @{$binTLen[$i]}, $aln[8];
                        #push @{$binLibs[$i]}, $libN;
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
		my $nFrags = (defined $binTLen[$i] ? @{$binTLen[$i]} : 0);
        my ($tLens, $libs);
		if ($nFrags == 0) { # bin has no data (R will ignore)
			($tLens, $libs) = ('NA', 'NA');
		} elsif ($nFrags <= $maxMlePairs) { # typical bin, use as is
			$tLens = join(",", @{$binTLen[$i]});
			#$libs  = join(",", @{$binLibs[$i]});
		} else { # downsample bins with too many pairs and proceed with KS/MLE
			my @is = (shuffle(0..$nFrags-1))[1..$maxMlePairs];
			$tLens = join(",", @{$binTLen[$i]}[@is]);
			#$libs  = join(",", @{$binLibs[$i]}[@is]);
            $nFrags = $maxMlePairs;
		}
        $prevChrom eq 'chrM' or
            #print $tLenH join("\t", $prevChrom, $binI, $nFrags, $tLens, $libs), "\n";
            print $tLenH join("\t", $prevChrom, $binI, $nFrags, $tLens, 'NA'), "\n";
			###################
			#$binI == 143746 and print join("\t", $prevChrom, $binI, $nFrags, $tLens, 'NA'), "\n";
    }
    splice(@binTLen, 0, $maxBin - $bin0 + 1);
	#splice(@binLibs, 0, $maxBin - $bin0 + 1);
    splice(@enc,     0, $maxBin - $bin0 + 1);
    $bin0 = $maxBin + 1;
}
sub finishChrom {
	my ($chrom) = @_;
	print STDERR "  finishing $prevChrom\n";	
	commitBins($rghtBinI);
	#($bin0, $leftBinI, $rghtBinI, @binTLen, @binLibs, @enc) = (1);
    ($bin0, $leftBinI, $rghtBinI, @binTLen, @enc) = (1);
	$prevChrom = $chrom;
}

# run MLE via R and commit the results
sub runMle {
	my ($tLenFile) = @_;
    $ENV{sample} = $sample;
	$ENV{libraries}     = join(",", sort {$libNs{$a} <=> $libNs{$b}} keys %libNs);
    $ENV{libraryNs}     = join(",", sort {$a <=> $b} values %libNs);    
    $ENV{pdfDir}        = $pdfDir;
    $ENV{tLenFile}      = $tLenFile;
    $ENV{minCnvSize}    = $minCnvSize;
	$ENV{maxCnvSize}    = $maxCnvSize;
    $ENV{cnvSizeStep}   = $cnvSizeStep;
    $ENV{min_mle_pairs} = $minMlePairs;
	$ENV{min_cnv_pairs} = $minCnvPairs;
    open(my $rH, "-|", "Rscript $libDir/mle.R") or die "could not open mle.R stream: $!\n";
    while (my $line = <$rH>) {
        chomp $line;
		my @f = split("\t", $line);
        #my ($chrom, $binI, $nFrags,
        #    $deltaCDF_ref, $deltaCDF_alt, $cnvSize, $cnvFrac, $cnvType, $cnvType2)
        my $bin = join("\t",
			$f[0], $f[1]*$genBinSize-1, $f[1]*$genBinSize,
			$sample, $f[2], '.',
			@f[3..8]
		)."\n";
        if ($region) { $f[2] and print STDOUT $bin } else { print $binH $bin }
		if ($prevBin and
			($$prevBin[0] ne $f[0] or
		     $$prevBin[8] ne $f[8])) { # merge dup and ins into single runs
			commitCNV();
		}
		$prevBin = \@f;
		($f[8] eq 'gain' or $f[8] eq 'loss') and push @cnvBins, \@f;
    }
	commitCNV();
    closeHandles($rH);
}
sub commitCNV {
	@cnvBins or return;
    unless($region){
        my $nBins = @cnvBins;
        if($nBins >= $minCnvBins){
            my $midI = int($nBins / 2);
            print $cnvH join("\t",
                $cnvBins[0][0], $cnvBins[0][1]*$genBinSize-1, $cnvBins[$#cnvBins][1]*$genBinSize,
                $sample, $nBins, '.',
                @{$cnvBins[$midI]}[3..8]
            ), "\n";			
        }
    }
	@cnvBins = ();	
}

#1034x1112p2
#[wilsonte@mengf-n6 bam4]$ echo $DUP2
#chr13:23392072-23419642
#[wilsonte@mengf-n6 bam4]$ echo $DEL
#chr11:88884000-88903000
#[wilsonte@mengf-n6 bam4]$ echo $REG
#chr11:88960000-88973000
#[wilsonte@mengf-n6 bam4]$

#OLD my ($chrom, $binI, $nFrags,
#    $nlog_ks_p_ref, nlog_ks_p_alt,
#    $dupSize, $nDupFrag, $delInsSize, $nDelInsFrag
#    $binType)
            
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
