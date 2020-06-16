use strict;
use warnings;
use List::Util qw(shuffle);

# define variables
use vars qw($version $utility $error $libDir $slurp
            $TRUE $FALSE
            $tmpH $inH $outH
            $inputDir $tmpDir $outputDir $plotDir);
require "$libDir/exclude_subs.pl";
my (@globs, $sample,
    $genomeFasta, $snpFile, $excludeFile, $maxMem,
    $minQual);
my $i = 0;
my %c = map { $_ => $i++ } qw(
    chrom start end name score strand
    refBase altBase    
    chrom_ start_ end_ readBases
);
my (%counts, %fracs);

# manage options
sub setOptions_genotype {
    setOptionValue(\$sample,     'sample');
    setOptionValue(\$genomeFasta,'genome-fasta');
    setOptionValue(\$snpFile,    'snp-file');
    setOptionValue(\$excludeFile,'exclude-file');  
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);
    setOptionValue(\$outputDir,  'output-dir');    
    setOptionValue(\$minQual,    'min-qual',        5);   
    #-------------------------------------
    -e $genomeFasta or die "$error: file not found: $genomeFasta\n";
    -e $snpFile or die "$error: file not found: $snpFile\n";
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    $minQual =~ m|\D| and die "$error: min-qual must be an integer MAPQ\n";    
}

# main execution block
sub svtools_genotype {
    (@globs) = @_;
    $globs[0] or die "$error: svtools genotype requires one or more coordinate-sorted bam files as input\n";

    # initialize
    setOptions_genotype();    
    print STDERR "$utility genotype: " . getTime(), "\n";
    initializeExclude($excludeFile);

    # prepare the mpileup version of the SNP file
    my $tmpFile = getTmpFile("snps", $sample, "txt");
    system("gunzip -c $snpFile | awk '{print \$1\"\\t\"\$3}' | $slurp -s 50M -o $tmpFile");
    
    # run mpileup against list of SNP positions
    my $inStream =
        "samtools merge -u - @globs | ".
        "samtools view -q $minQual -h - | ". # purge read groups, i.e. merge into a single sample
        "grep -v '\@RG' | ".
        "sed 's/\\tRG:Z:\\S*//' | ".
        "samtools view -Su - | ".
        "samtools rmdup - - 2>/dev/null | ".
        "samtools mpileup -f $genomeFasta -l $tmpFile - | ".
        "awk 'BEGIN{OFS=\"\\t\"}{print \$1,\$2-1,\$2,\$5,\$4,\".\"}' | ".
        "bedtools intersect -loj -a $snpFile -b - ";
    openInputStream($inStream, \$inH);

    # count the number of reads with each SNP value by position
	my $gntFile = getOutFile("snps", $sample, "bgz");
    openOutputStream($gntFile, \$outH, undef, undef,
            "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        my %b = (ref=>0, alt=>0);        
        if ($f[$c{chrom_}] ne '.') {
            my @b = split("", uc($f[$c{readBases}]));
            foreach my $i(0..$#b){
                $b[$i] or next;
                if ($b[$i] eq '^') { # e.g. ".^], 2"
                    $b[$i+1] = "";
                    next;
                }
                if ($b[$i] eq '+' or $b[$i] eq '-') { # e.g. ",..+4TCTC. 4" or "g.-2TG  2"
                    my $j = $i+1;
                    my $N = "";
                    while ($b[$j] =~ m/\d/) { $N .= $b[$j]; $j++ }
                    foreach my $j($i+1..$i+length($N)+$N) { $b[$j] = "" }
                    next;
                }            
                $b[$i] eq uc($f[$c{altBase}]) and $b{alt}++;
                ($b[$i] eq '.' or $b[$i] eq ',') and $b{ref}++;    
            }            
        }
        my $nReads = $b{ref} + $b{alt};
        print $outH join("\t",
            @f[$c{chrom}, $c{start}, $c{end}, $c{name}],
            $nReads, '.',
            @f[$c{refBase}, $c{altBase}],
            $b{ref}, $b{alt}, $f[$c{readBases}]
        ), "\n";
        $counts{$nReads}++;
        if ($nReads) {
            my $frac = int(($b{alt} / $nReads) / 0.05 + 0.5) * 0.05;
            $fracs{$frac}++
        }
    }
	closeHandles($inH, $outH);
    
    # report tallies
    printGenotypeSummary(\%counts, 'coverage',    qw(coverage nSNPs));
    printGenotypeSummary(\%fracs,  'frequencies', qw(fracAlt  nSNPs));
    
    # finish up
    unlink $tmpFile;
    system("tabix -p bed $gntFile");
}

sub printGenotypeSummary {
    my ($hash, $name, @lbls) = @_;
	my $outFile = getOutFile($name, $sample, "txt");
    openOutputStream($outFile, \$outH);    
    print $outH join("\t", @lbls), "\n";
    foreach my $val(sort {$a <=> $b} keys %$hash){
        print $outH join("\t", $val, $$hash{$val}), "\n";
    }
    closeHandles($outH);
}

#chr19   3078938 3078939 .       0       .       A       T       chr19   3078938 3078939 .       1       .
#chr19   3078946 3078947 .       0       .       C       T       chr19   3078946 3078947 .,      2       .
#chr19   3078961 3078962 .       0       .       G       T       chr19   3078961 3078962 ,,      2       .
#chr19   3078965 3078966 .       0       .       C       T       chr19   3078965 3078966 ,,      2       .
#chr19   3078985 3078986 .       0       .       G       T       chr19   3078985 3078986 ,,,     3       .
#chr19   3079044 3079045 .       0       .       G       C       chr19   3079044 3079045 .c.,,   5       .
#chr19   3079046 3079047 .       0       .       G       T       chr19   3079046 3079047 t.,,    4       .
#chr19   3094431 3094432 .       0       .       G       A       .       -1      -1      .       -1      .
#chr19   3094466 3094467 .       0       .       A       G       .       -1      -1      .       -1      .
#chr19   3099034 3099035 .       0       .       T       C       chr19   3099034 3099035 .,      2       .
#chr19   3099108 3099109 rs49045175      0       .       G       A       chr19   3099108 3099109 ,       1       .
#chr19   3099135 3099136 .       0       .       G       A       chr19   3099135 3099136 ,       1       .
#chr19   3099137 3099138 .       0       .       A       T       chr19   3099137 3099138 ,       1       .
#chr19   3099144 3099145 rs583749385     0       .       C       G       chr19   3099144 3099145 ,       1       .

#a dot stands for a match to the reference base on the forward strand,
#a comma for a match on the reverse strand,
#a '>' or '<' for a reference skip,
#`ACGTN' for a mismatch on the forward strand and
#`acgtn' for a mismatch on the reverse strand.

#A pattern `\\+[0-9]+[ACGTNacgtn]+'
#indicates there is an insertion between this reference position and the next reference position.
#The length of the insertion is given by the integer in the pattern, followed by the inserted sequence.
#Similarly, a pattern `-[0-9]+[ACGTNacgtn]+' represents a deletion from the reference.
#The deleted bases will be presented as `*' in the following lines.

#Also at the read base column, a symbol `^' marks the start of a read.
#The ASCII of the character following `^' minus 33 gives the mapping quality.
#A symbol `$' marks the end of a read segment.

            
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
