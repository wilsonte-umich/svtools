use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $_gap $_split $_clip);
our %idxSubs = ($_gap   => \&index_gap,
                $_split => \&index_split,
                $_clip  => \&index_clip);
our ($idxH, $libKey);

#========================================================================
# index 'extract'
#------------------------------------------------------------------------
# convert find format to tabix/bed format for web browser, etc.
sub index_gap {
    my ($chr1, $str1, $chr2, $str2,
        $outPos1, $maxDst1, $minDst1,
        $outPos2, $maxDst2, $minDst2,
        $innPos1, $innSde1, $innSeq1, $seq1,
        $innPos2, $innSde2, $innSeq2, $seq2) = @_;
    index_span( $_gap,
                $chr1, $chr2,
                $outPos1, $outPos2,
                $str1, $str2);
    $innPos1 and index_clip($chr1, $outPos1, $innSde1, $innSeq1, $seq1);
    $innPos2 and index_clip($chr2, $outPos2, $innSde2, $innSeq2, $seq2);
}
sub index_split {
    my ($chr1, $str1, $chr2, $str2,
        $idxPos1, $maxDst1, $minDst1,
        $idxPos2, $maxDst2, $minDst2,
        $jxnPos1, $jxnSde1,
        $jxnPos2, $jxnSde2,
        $microH, $insSeq, $seq) = @_;
    index_span( $_split,
                $chr1, $chr2,
                $jxnPos1, $jxnPos2,
                $str1, $str2);
}
sub index_span {
    my ($type,
        $chr1, $chr2,
        $pos1, $pos2,
        $str1, $str2) = @_;
    if ($chr1 ne $chr2) {
        print $idxH join("\t", $chr1, $pos1 - 1, $pos1, "$libKey,trans,$chr2,$pos2", $type, "+"), "\n";
        print $idxH join("\t", $chr2, $pos2 - 1, $pos2, "$libKey,trans,$chr1,$pos1", $type, "+"), "\n";
    } elsif($str1 == $str2){
        if ($str1) {
            print $idxH join("\t", $chr1, $pos1 - 1, $pos2, "$libKey,invR", $type, "+"), "\n";
        } else {
            print $idxH join("\t", $chr1, $pos1 - 1, $pos2, "$libKey,invF", $type, "+"), "\n";
        }
    } else {
        if ($pos1 <= $pos2){
            print $idxH join("\t", $chr1, $pos1 - 1, $pos2, "$libKey,del", $type, "+"), "\n";
        } else {
            print $idxH join("\t", $chr1, $pos2 - 1, $pos1, "$libKey,dup", $type, "+"), "\n";
        }
    }   
}
sub index_clip {
    my ($chr, $outPos, $outSde,
        $outSeq, $seq) = @_;
    print $idxH join("\t", $chr, $outPos - 1, $outPos, "$libKey,$outSde,$outSeq", $_clip, "+"), "\n";
}
#========================================================================

1;
