use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            $_gap $_split $_clip);
my (@globs, $region);
my ($chrom, $minPos, $maxPos);
my ($genomeFasta, $padding, $lineLen, $clippedOnly, $minQual);
my $lftClpRex = qr/^(\d+)(S|H)/;      
my $rgtClpRex = qr/(\d+)(S|H)$/;

# manage options
sub setOptions_pull {
    setOptionValue(\$genomeFasta,   'genome-fasta');
    setOptionValue(\$padding,       'padding',       0);
    setOptionValue(\$lineLen,       'line-length',   100);
    setOptionValue(\$clippedOnly,   'clipped-only');
    setOptionValue(\$minQual,       'min-qual',     0);
    $maxPos or $maxPos = $minPos;
    $minPos -= $padding;
    $maxPos += $padding;
    $region = "$chrom:$minPos-$maxPos";
}

# main execution block
sub svtools_pull {

    # parse the command line arguments
    $region = pop @_;
    $region =~ s/,//g;
    $region or die "$error: missing region as last argument\n";
    ($chrom, $minPos, $maxPos) = ($region =~ m/(\w+):(\d+)-(\d+)/);
    $chrom or ($chrom, $minPos) = ($region =~ m/(\w+):(\d+)/);
    $chrom or die "$error: invalid region as last argument: format chr:start[-end]\n";
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";

    # initialize
    setOptions_pull();    
    print STDERR "$utility pull: " . getTime(), "\n";

    # do the work
    pullReads();   
}

# coordinate the retrieval and display formatting of the reads in a region
sub pullReads {
    
    # pull the reads
    my $mrgStream = "samtools merge -u -R $region - @globs | samtools view -";
    open(my $inH, "-|", $mrgStream) or die_("could not open input stream: $!\n");
    my %alns;
    print "\n";
    print join("\t", qw(POS CIGAR TLEN SEQ)), "\n";
    while (my $line = <$inH>) {
        chomp $line;
        my @aln = split("\t", $line);
        $aln[1] & 2 and next; # a proper pair (only show anomalies)
        $aln[1] & 4 and next; # unmapped
        $aln[1] & 8 and next; # mate is unmapped
        $aln[4] < $minQual and next; # low quality read
        #$aln[6] eq '=' or next;
        my $pos; # the outermost position of the alignment IF all bases had mapped, including clips
        if($aln[1] & 0x010){
            $pos = getEnd($aln[3], $aln[5]) + ($aln[5] =~ m/$rgtClpRex/ ? $1 : 0);
        } else {
            $pos = $aln[3] - ($aln[5] =~ m/$lftClpRex/ ? $1 : 0);
        }
        my $key = join(":", $aln[2], $aln[1] & 0x010, $pos, $aln[1] & 0x040,
                            $aln[6], $aln[1] & 0x020, $aln[7]); # key to purge dups
        push @{$alns{$key}}, \@aln;
    }
    close $inH;
    print "\n";
    
    # parse the reads into sequence strings
    my @reads;
    foreach my $key(keys %alns){
        my $aln = (sort { $$b[4] <=> $$a[4] } @{$alns{$key}})[0]; # highest quality unique alignment
        my ($flag, $alnS, $cigar, $read) = @$aln[1, 3, 5, 9];
        (!$clippedOnly or $cigar =~ m/(S|H)/) or next;
        print join("\t", @$aln[3,5,8,9]), "\n";        
        my @read = split("", $read);
        my (@lftClip, @bases, @rgtClp, $bases);
        if ($cigar =~ m/^(\d+)(S|H)/) {
            $cigar = substr($cigar, length($1) + 1);
            @lftClip = $2 eq 'S' ? splice(@read, 0, $1) : "H" x $1;
        }
        if ($cigar =~ m/(\d+)(S|H)$/) {
            $cigar = substr($cigar, 0, -(length($1) + 1));
            @rgtClp = $2 eq 'S' ? splice(@read, -$1) : "H" x $1;
        }
        while ($cigar =~ (m/(\d+)(\w)/g)) {
            foreach my $i(0..$1-1){
                my $base = shift @read;
                $2 eq "I" or push @bases, $base;
            } 
        }
        my $lftClp = lc(join("", "", @lftClip));
        my $rgtClp = lc(join("", "", @rgtClp));
        my $alnE = getEnd($alnS, $cigar);
        my $plpS = $alnS - length($lftClp);
        my $plpE = $alnE + length($rgtClp);
        if($flag & 0x010){
            $plpS--;
            $bases = join("", "<", $lftClp, uc(join("", @bases)), $rgtClp);
        } else {
            $plpE++;
            $bases = join("", $lftClp, uc(join("", @bases)), $rgtClp, ">");
        }        
        $plpS < $minPos and next;
        $plpE > $maxPos and next;        
        push @reads, [$plpS, $plpE, join("", $bases)];
    }
    
    # report nothing to do
    unless(@reads){
        print STDERR "no anomalous reads to report in region\n";
        exit;
    }
    
    # get the genome sequence over which reads were returned
    my $refSeq = qx/samtools faidx $genomeFasta $region | tail -n+2/;
    chomp $refSeq;
    $refSeq =~ s/\s//g;
    my $plpLen = length($refSeq);
    my @lines = $refSeq;    
    
    # assemble the pileup lines
    my ($maxLine, %maxEnd) = (1);
    foreach my $read(sort {$$a[0] <=> $$b[0]} @reads){
        my ($plpS, $plpE, $aln) = @$read;
        my $line = 1;
        while ($plpS <= ($maxEnd{$line} || 0)) { $line++ }
        $maxEnd{$line} = $plpE;
        $lines[$line] or $lines[$line] = '-' x $plpLen;
        my $mapS = $plpS - $minPos;
        $lines[$line] = substr($lines[$line], 0, $mapS).
                        $aln.
                        substr($lines[$line], $mapS + length($aln));
    }
    
    # print the pileup
    my $lineStart = $minPos;
    print "\n";
    while ($lines[0]) {
        print "$lineStart\n";
        print "|\n";
        $lineStart += $lineLen;
        foreach my $i(0..$#lines){
            my $len = length($lines[$i]) < $lineLen ? length($lines[$i]) : $lineLen;
            my $line = substr($lines[$i], 0, $len);
            $lines[$i] = substr($lines[$i], $len);
            print "$line\n";
        }
        print "\n";
    }
    
    # convenience utility
    sub die_ {
        my ($msg) = @_;
        print "\n$error: $msg\n\n";
        exit;
    }    
}

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
