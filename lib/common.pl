use strict;
use warnings;
use Scalar::Util qw(looks_like_number);

#========================================================================
# 'common.pl' has miscellaneous subs common to many commands and utilities
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($utility $command @args  
            $slurp $genomeFasta);
our ($_gap, $_split, $_clip) = 0..2;
our $pTLenBin = 10;
our ($TRUE, $FALSE) = (1, undef);
our ($inH, $tmpH, $outH, $rptH);
our ($inputFile, $region, $ouputFile,
	 $inputDir, $tmpDir, $outputDir, $plotDir);
#========================================================================

#========================================================================
# sequence manipulation
#------------------------------------------------------------------------
sub revComp {
    my $seq = reverse(uc($_[0]));
    $seq =~ tr|ACGT|TGCA|;
    return $seq;
}
sub getRegSeqs {
	@_ or return [];
	open(my $faH, "-|", "samtools faidx $genomeFasta ".join(" ", @_))
		or die "could not open samtools faidx stream: $!\n";
	my $i = -1;
	my @regSeqs;
	while (my $line = <$faH>) {
		chomp $line;
		if ($line =~ m|^>|) { $i++ } else { $regSeqs[$i] .= uc($line) }
	}
	close $faH;
	return \@regSeqs;
}
#========================================================================

#========================================================================
# miscellaneous
#------------------------------------------------------------------------
sub reportCount {
    my ($n, $desc) = @_;
    print STDERR join("\t", sprintf("$utility $command: %13s", commify($n)), $desc), "\n";
}
sub commify {
    local $_  = shift;
    1 while s/^(-?\d+)(\d{3})/$1,$2/;
    return $_;
}
sub min {
    my ($v1, $v2) = @_;
    $v1 <= $v2 ? $v1 : $v2;
}
sub max {
    my ($v1, $v2) = @_;
    $v1 >= $v2 ? $v1 : $v2;
}
sub median {
    my (@data) = sort {$a <=> $b} @_;
    my $i = @data / 2;
    my $upper = $data[$i];
    @data % 2 and return $upper;
    my $lower = $data[$i - 1];    
    return($lower + ($upper - $lower) / 2);
}
sub mean{
    my (@data) = @_;
    @data == 0 and die "no values sent to mean\n";
    my $sum = 0;
    foreach (@data) { $sum += $_ }
    return $sum / @data;
}
sub stdev{
    my (@data) = @_;
    @data == 1 and return 0;
    my $mean = mean(@data);
    my $sqSum = 0;
    foreach(@data) { $sqSum += ($mean-$_) ** 2 }
    return ($mean, ($sqSum / (@data-1)) ** 0.5);
}
sub roundCount {
    my ($count, $scalar) = @_;
    $scalar or $scalar = 1000;
    int($count * $scalar + 0.5) / $scalar;
}
sub roundCount2 {
    my ($count) = @_;
    int($count * 100 + 0.5) / 100;
}
sub encode_json{
    my ($hash) = @_;
    my @json = ();
    foreach my $key(keys %$hash){
        my $quote = looks_like_number($$hash{$key}) ? "" : '"';
        push @json, "\"$key\" : $quote$$hash{$key}$quote"; 
    }
    return '\'{'.join(",",@json).'}\'';
}
#========================================================================

#========================================================================
# SAM handling functions
#------------------------------------------------------------------------  
sub getStats { # determine read statistics from first 100K SAM reads
	my ($inH, $repeat, $rptH, $statsLines, $nStatPairs, $minQual) = @_;	
    my ($inHeader, $n) = (1, 0);
    my ($minTLen, $maxTLen, $maxRLen) = (1e9, 0, 0);
    my ($libType, @tLens, @rLens, %tLens, %rLens, %pTLen, %chromSizes);
    while (my $line = <$inH>) {
        CHECKHEADER: if ($inHeader) {
            if ($line =~ m/^\@/) {
                $repeat and print $rptH $line; # only repeat header lines here
                $line =~ m/^\@SQ\tSN:(\S+)\tLN:(\d+)/ and $chromSizes{$1} = $2;  
            } else {
                $inHeader = undef;
                goto CHECKHEADER;   				
			}
        } else {
            $$statsLines .= $line;
            my @aln = split("\t", $line, 11);
            if(($aln[1] & 0x952) == 0x042 and # proper, first read, forward, not secondary or supplemental
                $aln[5] !~ m/S/ and           # unclipped
				$aln[4] >= $minQual and       # of sufficient uniqueness
				!checkExclude($aln[2], $aln[3] - 1, $aln[3] - 1 + $aln[8])){ # not excluded 
                $libType or $libType = ($aln[7] - $aln[3] > 0) ? "FR" : "RF";
				$minTLen <= $aln[8] or $minTLen = $aln[8];
                $maxTLen >= $aln[8] or $maxTLen = $aln[8];
                push @tLens, $aln[8];
				$tLens{$aln[8]}++;
                $pTLen{int($aln[8]/$pTLenBin + 0.5) * $pTLenBin}++;
                my $rLen = length($aln[9]);
                $maxRLen >= $rLen or $maxRLen = $rLen;
                push @rLens, $rLen;
				$rLens{$rLen}++;
                $n++; 
                last if($n >= $nStatPairs);
            }
        }
    }
    my $medTLen = (sort {$a <=> $b} @tLens)[int($#tLens/2)];
    my $medRLen = (sort {$a <=> $b} @rLens)[int($#rLens/2)];
    my $minPTLen = int($minTLen/$pTLenBin + 0.5) * $pTLenBin;
    my $maxPTLen = int($maxTLen/$pTLenBin + 0.5) * $pTLenBin;    
    my ($pShorter, $pLonger, $pTLen) = (0, 1, '');
    for (my $tLen = $minPTLen; $tLen <= $maxPTLen; $tLen += $pTLenBin){
        my $pTL = ($pTLen{$tLen} || 0) / $n;        
        $pShorter += $pTL;        
        $pTLen .= "$tLen=".min($pShorter, $pLonger).";";
        $pLonger -= $pTL;        
    }
	return ($libType, $minTLen, $maxTLen, $medTLen, $maxRLen, $medRLen,
            \%tLens, \%rLens, $pTLen, \%chromSizes);
}
sub getEnd { # rightmost mapped read position in reference genome
    my ($start, $cigar) = @_;
    $cigar =~ s/\d+(S|H)//g;
    my $end = $start - 1;
    while ($cigar =~ (m/(\d+)(\w)/g)) {
        $2 eq "I" or $end += $1;    
    }
    return $end; 
}
#========================================================================

#========================================================================
# anomalous pair criteria
#------------------------------------------------------------------------
sub setMinTLen {
    my ($caller, $libKey, $minTLenOpt, $minTLen, $tags) = @_;
    my ($minDel, $minDup, $minInv) = split("/", $minTLenOpt);
    %$minTLen = (
        del => int($$tags{maxTLen} + $minDel * ($$tags{maxTLen} - $$tags{medTLen})),
        dup => int($$tags{maxTLen} + $minDup * ($$tags{maxTLen} - $$tags{medTLen})),
        inv => int($$tags{maxTLen} + $minInv * ($$tags{maxTLen} - $$tags{medTLen})) 
    );
    print STDERR "$utility $caller: library = $libKey\n";
    print STDERR "$utility $caller: minDel = $$minTLen{del} bp\n";
    print STDERR "$utility $caller: minDup = $$minTLen{dup} bp\n";
    print STDERR "$utility $caller: minInv = $$minTLen{inv} bp\n"; 
}
sub checkTLen {
    my ($pairType, $rest, $minTLen) = @_;
    if ($pairType == $_gap) { # enforce minimum TLEN by pair-type and library
        my ($chr1, $str1, $chr2, $str2, $pos1, $xx, $yy, $pos2) = split("\t", $rest, 9);
        if ($chr1 eq $chr2) {
            my $tLen = abs($pos2 - $pos1);
            if ($str1 == $str2) { # inversion
                $tLen > $$minTLen{inv} or return;
            } elsif(($str2 and $pos2 > $pos1) or
                    ($str1 and $pos1 > $pos2)) { # deletion
                $tLen > $$minTLen{del} or return;
            } else { # duplication
                $tLen > $$minTLen{dup} or return;
            }
        }
    }
    return 1;
}
#========================================================================

1;
