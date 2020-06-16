use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir);
require "$libDir/exclude_subs.pl";
my (@globs, $sample, $library, $bam,
    $repeat, $nStatPairs, $excludeFile, $minQual);
my ($stepsAfter, $stepsBefore);
my ($libKey, $statsLines, %tLens, $p, %chromStats);
my ($libType, $minTLen, $maxTLen, $maxRLen, $medTLen, $medRLen,
    $tLens, $rLens, $pTLen, $chromSizes) = (1e9, 0, 0);
my ($pos0, $nFrags, $sLen, @breaks) = (0, 0, 0); 

# manage options
sub setOptions_map {
    setOptionValue(\$sample,    'sample');
    setOptionValue(\$library,   'library');       
    setOptionValue(\$bam,       'bam');
    setOptionValue(\$tmpDir,    'tmp-dir',     '/tmp');
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$repeat,    'repeat');    
    setOptionValue(\$nStatPairs,'n-stat-pairs', 1e5);
    setOptionValue(\$excludeFile,'exclude-file');
    setOptionValue(\$minQual,   'min-qual',     5);
    #-------------------------------------
    $libKey = "$sample.$library";
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $nStatPairs  =~ m|\D| and die "$error: option --n-stat-pairs must be an integer number\n";
    $minQual =~ m|\D| and die "$error: min-qual be an integer MAPQ\n";
    if ($bam) {
        $stepsAfter .= "samtools view -h -";
        $repeat and $stepsBefore = "samtools view -hSb -";
    }
}

# main execution block
sub svtools_map {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";
    $globs[1] and die "$error: only one input stream or file allowed for command 'map'; ".
                      "use 'samtools merge ... | svtools map [options] -' to merge multiple bam files\n";
    $globs[0] eq '-' or $globs[0] eq 'stdin' or -f $globs[0] or die "$error: unknown file: $globs[0]\n";

    # initialize
    setOptions_map();    
    print STDERR "$utility map: " . getTime(), "\n";
    initializeExclude($excludeFile);

    # open output streams
    my $tmpFile = getTmpFile("breaks", $libKey, "gz");
    openOutputStream($tmpFile, \$tmpH, $TRUE);
    $repeat and openOutputStream(undef, \$rptH, undef, undef, $stepsBefore);

    # do the work
    processInputFiles(\$inH, \@globs, $TRUE, $stepsAfter, $FALSE, \&handle_map_file);
    closeHandles($tmpH, $rptH);
    
    # commit the final file
    openInputStream($tmpFile, \$inH, $TRUE);
    my $outFile = getOutFile("breaks", $libKey, "gz");
    openOutputStream($outFile, \$outH, $TRUE);
    my (@chromSizes, @chromStats);
    foreach my $chrom(sort keys %$chromSizes){
        $chromStats{$chrom} or next;
        push @chromSizes, "$chrom,$$chromSizes{$chrom}";
        push @chromStats, "$chrom,$chromStats{$chrom}";
        print STDERR "$utility map: $chrom $chromStats{$chrom}\n";
    }
    my $chromSizes = join(";", @chromSizes);
    my $chromStats = join(";", @chromStats);
    print $outH join("\t", "#application:svtools",
                            "version:$version",
                            "file_type:map",
                            "sample:$sample",
                            "library:$library",
                            "libType:$libType",
                            "maxTLen:$maxTLen",
                            "medTLen:$medTLen",
                            "medRLen:$medRLen",
                            "nStatPairs:$nStatPairs",
                            "chromSizes:$chromSizes",
                            "chromStats:$chromStats"), "\n";
    while (<$inH>) { print $outH $_ }
    closeHandles($inH, $outH);
    unlink $tmpFile;

    # report counts
    reportCount($libType,    "$libKey\tlibrary type");
    reportCount($medRLen,    "$libKey\tmedian read length");
    reportCount($maxRLen,    "$libKey\tlargest read length");
    reportCount($medTLen,    "$libKey\tmedian expected insert size");
    reportCount($minTLen,    "$libKey\tsmallest expected insert size");
    reportCount($maxTLen,    "$libKey\tlargest expected insert size");
}

# process each input file
my $fileN = 1;
sub handle_map_file {
    my ($file) = @_;
    print STDERR "$utility map: $file\n";
    
    # get pair statistics
    if ($fileN == 1) {
        print STDERR "$utility map: getting read pair statistics\n";
        ($libType, $minTLen, $maxTLen, $medTLen, $maxRLen, $medRLen, $tLens, $rLens, $pTLen, $chromSizes) =
            getStats($inH, $repeat, $rptH, \$statsLines, $nStatPairs, $minQual);
        print STDERR "$utility map: mapping expected fragment coverage\n";
        open(my $statsH, "<", \$statsLines) or die "$error: could not open handle to stats lines: $!\n";
        run_map($statsH);
        close $statsH;
        $statsLines = ""; # release memory  
    }

    # run break mapping to tmp file
    eof($inH) or run_map($inH);
    closeHandles($inH);
    processPosGroup();
    finishChrom(); # must be here so stats pairs don't finish chrom prematurely
    $fileN++;
}

# parse the break points    
sub run_map {
    my ($inH_) = @_;
    while(my $line = <$inH_>){
        $line =~ m/^\@/ and next;
        $repeat and print $rptH $line;
        my @f = split("\t", $line, 10);
        ($f[1] & 0x912) == 0x002 or next; # forward read of a proper pair, not secondary or supplemental
        $f[8] > 0 or next;                # this is the leftmost read of a pair (i.e F/R)
        $f[4] >= $minQual or next;        # or sufficent MAPQ
        checkExclude($f[2], $f[3] - 1, $f[3] - 1 + $f[8]) and next; # have been instructed to ignore
        if ($p and ($f[2] ne $$p[2] or $f[3] != $$p[3])) {
            processPosGroup();
            $f[2] ne $$p[2] and finishChrom();
        }
        $tLens{$f[8]}++; # hash purges fragment dups based on length from each chrom/pos
        $p = \@f;
    }
}

sub finishChrom {
    printBreaks($#breaks);
    $chromStats{$$p[2]} = join(",",
        $nFrags,
        int($sLen / $nFrags + 0.5), # averageLen
        int($sLen / $$chromSizes{$$p[2]} * 10 + 0.5) / 10 # fragment coverage
    );
    ($pos0, $nFrags, $sLen, @breaks) = (0, 0, 0); 
}
sub processPosGroup {
    my @lens = keys %tLens; 
    $pos0 and printBreaks($#breaks < $$p[3]-$pos0-1 ? $#breaks : $$p[3]-$pos0-1);
    $pos0 = $$p[3]; 
    $breaks[0] += @lens; # increment the coverage plot up
    foreach my $len(@lens){
        $sLen += $len;
        $nFrags++;
        $breaks[$len]--; # increment the coverage plot down
    }
    %tLens = ();
}
sub printBreaks {
    my ($maxI) = @_;
    foreach my $i(0..$maxI){
        my $n = shift @breaks;
        $n or next;
        print $tmpH join("\t", $$p[2], $i + $pos0, $n), "\n"; # chrom, pos, increment
    }  
}

#Col	Field	Description
#1	QNAME	Query template/pair NAME
#2	FLAG	bitwise FLAG
#3	RNAME	Reference sequence NAME
#4	POS	1-based leftmost POSition/coordinate of clipped sequence
#5	MAPQ	MAPping Quality (Phred-scaled)
#6	CIAGR	extended CIGAR string
#7	MRNM	Mate Reference sequence NaMe (= if same as RNAME)
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
