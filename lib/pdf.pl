use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inputDir $outputDir $plotDir);
require "$libDir/exclude_subs.pl";
my (@globs, $sample, $excludeFile,
    $minQual, $maxDistTLen, $tLenBinSize, $nStatPairs);
my ($bamH, $distH);
my ($libraryN, %libs, %libNs, %distFiles) = (0);

# manage options
sub setOptions_pdf {
    setOptionValue(\$sample,     'sample');
    setOptionValue(\$excludeFile,'exclude-file');    
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$plotDir,    'plot-dir');    
    setOptionValue(\$minQual,    'min-qual',     5);
    setOptionValue(\$maxDistTLen,'max-tLen',     20000);
    setOptionValue(\$tLenBinSize,'tLen-bin-size',200);
    setOptionValue(\$nStatPairs, 'n-stat-pairs', 1e6);
    #-------------------------------------
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    $plotDir and (-d $plotDir or die "$error: directory not found: $plotDir\n");
    $minQual =~ m|\D| and die "$error: min-qual must be an integer MAPQ\n";
	$maxDistTLen =~ m|\D| and die "$error: max-tLen must be an integer number of bp\n";
    $tLenBinSize =~ m|\D| and die "$error: tLen-bin-size must be an integer number of bp\n";
}

# main execution block
sub svtools_pdf {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";

    # initialize
    setOptions_pdf();    
    print STDERR "$utility pdf: " . getTime(), "\n";
    initializeExclude($excludeFile);
    
    # determine pair statistics for all merged libraries of a single sample
	#	!! expects just one coordinate-sorted bam file per library, with library specified in RG tag
	# 	!! expects all files/libraries to come from just one source sample	    
    print STDERR "determining TLEN cumulative distribution frequencies\n";    
    if($globs[1]){
        my $mrgStream = "samtools merge -u - @globs | samtools view -";
        openInputStream($mrgStream, \$bamH);
        getPDF('merge of all input files', 'all');
    }
    
	# determine pair statistics for each input library for this sample
    processInputFiles(\$bamH, \@globs, $FALSE, "samtools view -", $FALSE, \&getPDF);
  
    # plot the reults
    my @libraries = sort { $libNs{$a} <=> $libNs{$b} } keys %libNs;
    my $libraries = join(",", @libraries);
    my $distFiles = join(",", map { $distFiles{$_} } @libraries);
    my $jpgFile   = getOutFile($sample, 'plot', "jpg");
    $ENV{sample}    = $sample;
    $ENV{libraries} = $libraries;
    $ENV{distFiles} = $distFiles;
    $ENV{jpgFile}   = $jpgFile;
    system("Rscript $libDir/pdf.R");
}

# calculate pair stats from the first n-stat-pairs pairs for a given library
sub getPDF {
	
	# assign a number to this file/library
    my ($file, $library) = @_;
	$libraryN++;
    $library and $libs{$libraryN} = $library;
	print STDERR $file, "\n";
	
	# extract sufficient pairs to determine FR (proper + deletion) TLEN distribution
	my ($n, %enc, %tLens) = (0);	
    while (my $line = <$bamH>) {
		my @aln = split("\t", $line, 12);		
		unless($libs{$libraryN}){
			$aln[11] =~ m/RG:Z:(\S+)/ or die "missing read group with library information in $file\n";
			$libs{$libraryN} = $1;
			%libNs = reverse %libs;
			print STDERR "  library = $libs{$libraryN}\n";
		}	
		if(($aln[1] & 0x930) == 0x020 and # forward read of FR pair, not secondary or supplemental
			$aln[6] eq '=' and            # on same chromosome
			$aln[8] >= -$maxDistTLen and  # within provided size range
            $aln[8] <=  $maxDistTLen and
			!$enc{"$aln[3],$aln[8]"} and  # not a duplicate of previously seen pair			
			$aln[4] >= $minQual and       # of sufficient uniqueness
			!checkExclude($aln[2], $aln[3] - 1, $aln[3]) ){ # not excluded	
				$enc{"$aln[3],$aln[8]"}++; # suppress duplicate pairs
				my $tLen = int($aln[8] / $tLenBinSize + 0.5) * $tLenBinSize;           
				$tLens{$tLen}++;
				$n++; 
				last if($n >= $nStatPairs);
		}
    }
    
    # report a bit of information
	my ($cdf, $tLen50, $tLen99) = (0);
	foreach my $tLen (sort {$a <=> $b} keys %tLens){
		$cdf += $tLens{$tLen} / $n;
        $tLen50 or ($cdf > 0.5  and $tLen50 = $tLen);
        $tLen99 or ($cdf > 0.99 and $tLen99 = $tLen);
		$cdf > 0.99 and last;
	}
	print STDERR "  TLEN at 50th percentile = $tLen50\n";	
	print STDERR "  TLEN at 99th percentile = $tLen99\n";	
    
	# establish PDF    
    my $distFile = getOutFile($sample, $libs{$libraryN}, "gz");
    $distFiles{$libs{$libraryN}} = $distFile;
    openOutputStream($distFile, \$distH, $TRUE);
    my $maxBin = int($maxDistTLen / $tLenBinSize + 0.5) * $tLenBinSize; 
    for (my $tLen=-$maxBin; $tLen<=$maxBin; $tLen+=$tLenBinSize){
        my $freq = ($tLens{$tLen} || 1) / $n;
        print $distH "$tLen\t$freq\n";
    }
    closeHandles($bamH, $distH);
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
