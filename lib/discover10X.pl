use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH
            $inputDir $tmpDir $outputDir $plotDir);
#require "$libDir/exclude_subs.pl";
my (@globs, $sample, $cellN, $minQual);
my ($weightPerCell, $chromWeight);
my %allowedChroms = map { ($_ => 1, "chr$_" => 1) } (1..100, "X", "Y");
map { $allowedChroms{$_} => 1} keys %allowedChroms;
my (%chroms, @chroms, %counts, %trueCells);
my $nTrueCells = 0;
my $f = 0x001 + # the read is paired in sequencing
        0x002;  # the read is mapped in a proper pair
my $F = 0x004 + # the query sequence itself is unmapped
        0x008 + # the mate is unmapped
        0x010 + # strand of the query (1 for reverse)
        0x100 + # the alignment is not primary
        0x200 + # the read fails platform/vendor quality checks
        0x400 + # the read is either a PCR or an optical duplicate
        0x800;  # supplementary alignment

# manage options
sub setOptions_discover10X {
    setOptionValue(\$cellN,     'expected-cells');
    setOptionValue(\$minQual,   'min-qual',     30);    
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$sample,    'sample');
    #-------------------------------------
    $cellN =~ m|\D| and die "$error: expected-cells must be an integer\n";
    $minQual =~ m|\D| and die "$error: min-qual must be an integer\n";
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
    $weightPerCell = 100;
    $chromWeight = $cellN * $weightPerCell;
}

# main execution block
sub svtools_discover10X {
    (@globs) = @_;
    my $usage = "discover10X only accepts the name of a 10X Cell Ranger DNA BAM file as input";
    $globs[0] or die "$error: no input file specified: $usage\n";
    ($globs[0] eq '-' or $globs[0] eq 'stdin') and die "$error: $usage\n";

    # initialize
    setOptions_discover10X();    
    print STDERR "$utility discover10X: " . getTime(), "\n";
    #initializeExclude($excludeFile);
    
    # do the work
    my $outFile = getOutFile("cells", "$sample", "txt");
    openOutputStream($outFile, \$outH);
    loadChromosomes10X();
    run_discover10X();
    print_discover10X();
    closeHandles($outH);
    
    # report counts
    reportCount(scalar(keys %counts),  "cell barcodes encountered in ".scalar(@chroms)." test regions");
    reportCount($nTrueCells,           "cells accepted as true cells");
}

# discover the chromosomes in use
sub loadChromosomes10X {
    print STDERR "$utility discover10X: loading chromosome information\n";
    openBamHeader($globs[0], \$inH);
    while (my $line = <$inH>){
        chomp $line;
        if($line =~ m/^\@SQ/){
            my ($chrom, $length) = ($line =~ m/SN:(\S+)\s+LN:(\d+)/) or die "invalid header SQ line: $line\n";
            $allowedChroms{$chrom} or next;            
            $chroms{$chrom} = $length;
        }
    }
    closeHandles($inH);
    @chroms = sort keys %chroms;
}

# count limited reads by barcode over all allowed chromosomes
sub run_discover10X {
    print STDERR "$utility discover10X: counting reads for barcodes\n";
    foreach my $chrom(@chroms){
        print STDERR "\tchrom: $chrom\n";
        my $start = int(rand(int($chroms{$chrom} * 0.8)));
        openBamStream($globs[0], \$inH, "50M", $minQual, $f, $F, "$chrom:$start");
        my $nReads = 0;
        my %chromCounts;
        while(my $line = <$inH>){
            my @f = split("\t", $line, 10);
            $f[8] > 0 or next; # this is the leftmost read of a pair (i.e F/R)
            $f[9] =~ m/CB:Z:(\w+)/ or next; # the 10X error corrected, validate cell barcode
            $counts{$1}{$chrom}++;
            $chromCounts{$1}++;
            $nReads++;
            $nReads >= $chromWeight and last;
        }
        closeHandles($inH);
        scoreTrueCells($chrom, \%chromCounts);
    }
}
sub scoreTrueCells { # currently the same heuristic as Cell Ranger DNA, but using a subset of reads per chrom
    my ($chrom, $chromCounts) = @_;
    my @CB = sort { $$chromCounts{$b} <=> $$chromCounts{$a} } keys %$chromCounts;
    my @N  = map { $$chromCounts{$_} } @CB;
    my $most = $N[0];
    my $threshold = $most * 0.1;
    my $i = 0;
    while($N[$i] >= $threshold){ $i++ }
    my $perc99 = $N[int($i * 0.01)];
    $threshold = $perc99 * 0.1;
    $i = 0;
    while($N[$i] >= $threshold){ $i++ }
    foreach my $j(0..$i){
        $trueCells{$CB[$j]}{$chrom} = 1;  
    }
}

# print a table of results by barcode
# table includes a row for every barcode accepted on at least one chromosome
# column "true_cell" == 1 if barcode was accepted on at least half of all test regions
sub print_discover10X {
    print $outH join("\t", qw(cell true_cell n_reads n_chroms), @chroms), "\n";
    foreach my $CB(keys %trueCells){
        my ($n_reads, $n_chroms, @chromResults) = (0, 0);
        foreach my $chrom(@chroms){
            my $trueCell = $trueCells{$CB}{$chrom} || 0;
            push @chromResults, $trueCell;
            $n_chroms += $trueCell;
            $n_reads += $counts{$CB}{$chrom} || 0;
        }
        my $trueCell = $n_chroms >= @chroms/2 ? 1 : 0;
        $nTrueCells += $trueCell;
        print $outH join("\t", $CB, $trueCell, $n_reads, $n_chroms, @chromResults), "\n";
    } 
}

1;

#...the following simple heuristic to identify cell barcodes:
#    we order the barcodes by the number of reads in ascending order
#    we select the barcodes that contain at least 1/10th as many reads as the last barcode. This gives us an initial set enriched for cell barcodes.
#    we calculate the 99-th percentile of the reads per barcode distribution restricted to the set of barcodes identified in the previous step
#    we define all barcodes with at least 1/10th as many reads as the 99-th percentile as cell barcodes

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

#@HD     VN:1.3  SO:coordinate
#@SQ     SN:1    LN:249250621
#@SQ     SN:2    LN:243199373
#@SQ     SN:3    LN:198022430
#@SQ     SN:4    LN:191154276
#@SQ     SN:5    LN:180915260
#@SQ     SN:6    LN:171115067
#@SQ     SN:7    LN:159138663
#@SQ     SN:8    LN:146364022
#@SQ     SN:9    LN:141213431
#@SQ     SN:10   LN:135534747
#@SQ     SN:11   LN:135006516
#@SQ     SN:12   LN:133851895
#@SQ     SN:13   LN:115169878
#@SQ     SN:14   LN:107349540
#@SQ     SN:15   LN:102531392
#@SQ     SN:16   LN:90354753
#@SQ     SN:17   LN:81195210
#@SQ     SN:18   LN:78077248
#@SQ     SN:19   LN:59128983
#@SQ     SN:20   LN:63025520
#@SQ     SN:21   LN:48129895
#@SQ     SN:22   LN:51304566
#@SQ     SN:X    LN:155270560
#@SQ     SN:Y    LN:59373566
#@SQ     SN:MT   LN:16569
#@SQ     SN:GL000207.1   LN:4262
#@SQ     SN:GL000226.1   LN:15008
#@SQ     SN:GL000229.1   LN:19913
#@SQ     SN:GL000231.1   LN:27386
#@SQ     SN:GL000210.1   LN:27682
#@SQ     SN:GL000239.1   LN:33824
#@SQ     SN:GL000235.1   LN:34474
#@SQ     SN:GL000201.1   LN:36148
#@SQ     SN:GL000247.1   LN:36422
#@SQ     SN:GL000245.1   LN:36651
#@SQ     SN:GL000197.1   LN:37175
#@SQ     SN:GL000203.1   LN:37498
#@SQ     SN:GL000246.1   LN:38154
#@SQ     SN:GL000249.1   LN:38502
#@SQ     SN:GL000196.1   LN:38914
#@SQ     SN:GL000248.1   LN:39786
#@SQ     SN:GL000244.1   LN:39929
#@SQ     SN:GL000238.1   LN:39939
#@SQ     SN:GL000202.1   LN:40103
#@SQ     SN:GL000234.1   LN:40531
#@SQ     SN:GL000232.1   LN:40652
#@SQ     SN:GL000206.1   LN:41001
#@SQ     SN:GL000240.1   LN:41933
#@SQ     SN:GL000236.1   LN:41934
#@SQ     SN:GL000241.1   LN:42152
#@SQ     SN:GL000243.1   LN:43341
#@SQ     SN:GL000242.1   LN:43523
#@SQ     SN:GL000230.1   LN:43691
#@SQ     SN:GL000237.1   LN:45867
#@SQ     SN:GL000233.1   LN:45941
#@SQ     SN:GL000204.1   LN:81310
#@SQ     SN:GL000198.1   LN:90085
#@SQ     SN:GL000208.1   LN:92689
#@SQ     SN:GL000191.1   LN:106433
#@SQ     SN:GL000227.1   LN:128374
#@SQ     SN:GL000228.1   LN:129120
#@SQ     SN:GL000214.1   LN:137718
#@SQ     SN:GL000221.1   LN:155397
#@SQ     SN:GL000209.1   LN:159169
#@SQ     SN:GL000218.1   LN:161147
#@SQ     SN:GL000220.1   LN:161802
#@SQ     SN:GL000213.1   LN:164239
#@SQ     SN:GL000211.1   LN:166566
#@SQ     SN:GL000199.1   LN:169874
#@SQ     SN:GL000217.1   LN:172149
#@SQ     SN:GL000216.1   LN:172294
#@SQ     SN:GL000215.1   LN:172545
#@SQ     SN:GL000205.1   LN:174588
#@SQ     SN:GL000219.1   LN:179198
#@SQ     SN:GL000224.1   LN:179693
#@SQ     SN:GL000223.1   LN:180455
#@SQ     SN:GL000195.1   LN:182896
#@SQ     SN:GL000212.1   LN:186858
#@SQ     SN:GL000222.1   LN:186861
#@SQ     SN:GL000200.1   LN:187035
#@SQ     SN:GL000193.1   LN:189789
#@SQ     SN:GL000194.1   LN:191469
#@SQ     SN:GL000225.1   LN:211173
#@SQ     SN:GL000192.1   LN:547496
#@SQ     SN:NC_007605    LN:171823
#@SQ     SN:hs37d5       LN:35477943
#@RG     ID:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1  SM:bj_mkn45_10pct       LB:MissingLibrary.1     PU:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1  PL:ILLUMINA
#@RG     ID:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:2  SM:bj_mkn45_10pct       LB:MissingLibrary.1     PU:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:2  PL:ILLUMINA
#@RG     ID:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:3  SM:bj_mkn45_10pct       LB:MissingLibrary.1     PU:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:3  PL:ILLUMINA
#@RG     ID:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:4  SM:bj_mkn45_10pct       LB:MissingLibrary.1     PU:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:4  PL:ILLUMINA
#@PG     PN:bwa  ID:bwa  VN:0.7.12-r1039 CL:bwa mem -p -t 4 -M -R @RG\tID:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1\tSM:bj_mkn45_10pct\tLB:MissingLibrary.1\tPU:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1\tPL:ILLUMINA /mnt/scratch2/abby/refdata-GRCh37-1.0.0/fasta/genome.fa /mnt/scratch2/abby/bj_mkn45_10pct/CNV_CALLER_SINGLECELL_CS/CNV_CALLER_SINGLECELL/_ALIGNER/TRIM_READS/fork0/chnk0-u1e792fc619/files/read1.fastq
#@PG     PN:longranger.attach_bcs        ID:attach_bcs   VN:1.0.0        PP:bwa
#@PG     PN:longranger.mark_duplicates   ID:mark_duplicates      VN:1.0.0        PP:attach_bcs
#@CO     10x_bam_to_fastq:R1(CR:CY,SEQ:QUAL)
#@CO     10x_bam_to_fastq:R2(SEQ:QUAL)
#@CO     10x_bam_to_fastq:I1(BC:QT)

#A00228:250:H3TVLDSXX:1:1146:8477:8469   1171    22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    NM:i:1  MD:Z:34G65      AS:i:95XS:i:90  XA:Z:14,+19789081,100M,2;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFFFFFF   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:1:1146:8160:8484   147     22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FFFFFFF,FFFFFFFF:FFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    NM:i:1  MD:Z:34G65      AS:i:95XS:i:90  XA:Z:14,+19789081,100M,2;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFFFFFF   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:1:2628:1479:16297  1171    22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCATGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FF,F:FFFF:FFFF:FF,,FF,FFFFFFFF,FFFFFFFFFFF,:::FFFFFFFF,,FF::FFF:FFFFFFFFFFFFFFFFFFFFFFFFF,,FFFF:FFFF    NM:i:2  MD:Z:34G7G57    AS:i:90XS:i:85  XA:Z:14,+19789081,100M,3;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFF,FFFFF,FF:FF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFF,F,F   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:2:1373:28854:34428 147     22      16053829        50      100M    =       16053669        -260   TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCGAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     :FFF:FFFFFF::FFF:FFF,FFF:FFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF    NM:i:0  MD:Z:100        AS:i:10XS:i:95  XA:Z:14,+19789081,100M,1;       CR:Z:CGGAGCTAGATCCGAG   CY:Z:FFF:FFFFFFFFFFFF
#CB:Z:CGGAGCTAGATCCGAG-1 BC:Z:CAATGGAG   QT:Z::::FF,,,   GP:i:698299001  MP:i:698298741  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:2
