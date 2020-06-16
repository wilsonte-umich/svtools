use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            $_gap $_split $_clip);
require "$libDir/exclude_subs.pl";
my (@globs, $sample, $library, $bam, $maxMem,
    $repeat, $excludeFile, $minQual);
my ($libKey, $stepsAfter, $stepsBefore);
my $lftClpRex = qr/^(\d+)(S|H)/;      
my $rgtClpRex = qr/(\d+)(S|H)$/;
my ($nReads, $nMasked) = (0) x 10;

# manage options
sub setOptions_mask {
    setOptionValue(\$sample,    'sample');
    setOptionValue(\$library,   'library');    
    setOptionValue(\$bam,       'bam');
    setOptionValue(\$tmpDir,    'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,    'max-mem',      1000000000);
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$repeat,    'repeat');    
    setOptionValue(\$excludeFile,'exclude-file');
    setOptionValue(\$minQual,   'min-qual',     5);
    #-------------------------------------
    $libKey = "$sample.$library";
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $plotDir and (-d $plotDir or die "$error: directory not found: $plotDir\n");
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";
    $minQual =~ m|\D| and die "$error: min-qual must be an integer MAPQ\n";
    if ($bam) {
        $stepsAfter .= "samtools view -h -";
        $repeat and $stepsBefore = "samtools view -hSb -";
    }
}

# main execution block
sub svtools_mask {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";

    # initialize
    setOptions_mask();    
    print STDERR "$utility mask: " . getTime(), "\n";
    initializeExclude($excludeFile);

    # open output streams
    my $outFile = getOutFile("knowns", $libKey, "bgz");
    my $bgZip = "sort -S $maxMem"."b -k1,1 -k2,2n | uniq | bgzip -c";    
    openOutputStream($outFile, \$outH, undef, undef, $bgZip);
    $repeat and openOutputStream(undef, \$rptH, undef, undef, $stepsBefore);
    
    # do the work
    processInputFiles(\$inH, \@globs, $TRUE, $stepsAfter, $FALSE, \&handle_mask_file);
    closeHandles($rptH, $outH, $inH);
    system("tabix -s 1 -b 2 -e 2 $outFile");

    # report counts
    reportCount($nReads,    "input reads");
    reportCount($nMasked,   "masked reads");
}

# process each input file
sub handle_mask_file {
    my ($file) = @_;
    print STDERR "$utility mask: $file\n";
    eof($inH) or run_mask($inH);
    closeHandles($inH);  
}

# create mask of all known ends from observed expected/germline/reference fragments
sub run_mask {
    my ($inH_) = @_;
    
    # get the alignment
    while (my $line = <$inH>) {
        $line =~ m/^\@/ and next;
        $repeat and print $rptH $line;
        $nReads++;
        
        # determine whether it should be considered
        my @aln = split("\t", $line, 7);
        $aln[4] >= $minQual or next; # a high quality read         
        $aln[1] & 0x002 or next;     # from a proper pair        
        $aln[1] & 0x900 and next;    # not secondary or supplemental
        checkExclude($aln[2], $aln[3] - 1, $aln[3]) and next; # not excluded
        $nMasked++;
        
        # save a line in the file for each encountered proper fragment end
        # such ends ought not to appear in true SV fragments
        
        # TODO: add UMIs to the read identification
        
        #my $pos; # the outermost position of the alignment IF all bases had mapped, including clips
        #if($aln[1] & 0x010){
        #    $pos = getEnd($aln[3], $aln[5]) + ($aln[5] =~ m/$rgtClpRex/ ? $1 : 0);
        #} else {
        #    $pos = $aln[3] - ($aln[5] =~ m/$lftClpRex/ ? $1 : 0);
        #}
        #print $outH join("\t", $aln[2], $pos,
        #                 $aln[1] & 0x010 ? "-" : "+",
        #                 $aln[1] & 0x040 ? 1 : 2), "\n"; # find is NOT currently using readN
        
        my $pos; # outPos as used by extract and find (not including clipped bases)
        if($aln[1] & 0x010){
            $pos = getEnd($aln[3], $aln[5]);
        } else {
            $pos = $aln[3];
        }
        print $outH join("\t", $aln[2], $pos, $aln[1] & 0x010), "\n"; # find is NOT currently using readN
    }
}

# 1111111^1111111>-------^----------------^---------
# -------^---------------^---------<222222^222222222

#Col	Field	Description
#1	QNAME	Query template/pair NAME
#2	FLAG	bitwise FLAG
#3	RNAME	Reference sequence NAME
#4	POS	1-based leftmost POSition/coordinate of clipped sequence
#5	MAPQ	MAPping Quality (Phred-scaled)
#6	CIAGR	extended CIGAR string
#7	MRNM	Mate Reference sequence NaMe
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
