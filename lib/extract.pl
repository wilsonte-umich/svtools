use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            $_gap $_split $_clip);
require "$libDir/exclude_subs.pl";
my (@globs, $sample, $library, $bam, $coordinate, $maxMem,
    $repeat, $nStatPairs, $excludeFile, $minClip, $minQual);
my ($gZip, $tmpBam, $stepsAfter, $stepsBefore);
my ($statsLines, $key, $libKey);
my ($libType, $minTLen, $maxTLen, $maxRLen, $medTLen, $medRLen, $tLens, $rLens, $pTLen) = (1e9, 0, 0);
my $prevName = "NA";
my @alns = ();
my $lftClpRex = qr/^(\d+)(S|H)/;      
my $rgtClpRex = qr/(\d+)(S|H)$/;
my %anom = (inv=>{}, del=>{}, dup=>{});
my ($nPairs, $nGaps, $nSplits, $nBadSplits, $nClips) = (0) x 100;

# manage options
sub setOptions_extract {
    setOptionValue(\$sample,    'sample');
    setOptionValue(\$library,   'library');    
    setOptionValue(\$bam,       'bam');
    setOptionValue(\$coordinate,'coordinate');
    setOptionValue(\$tmpDir,    'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,    'max-mem',      1000000000);
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$plotDir,   'plot-dir');
    setOptionValue(\$repeat,    'repeat');    
    setOptionValue(\$nStatPairs,'n-stat-pairs', 1e5);
    setOptionValue(\$excludeFile,'exclude-file');
    setOptionValue(\$minClip,   'min-clip',     8);
    setOptionValue(\$minQual,   'min-qual',     5);
    #-------------------------------------
    $libKey = "$sample.$library";
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $plotDir and (-d $plotDir or die "$error: directory not found: $plotDir\n");
    $nStatPairs  =~ m|\D| and die "$error: option --n-stat-pairs must be an integer number\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";
    $minQual =~ m|\D| and die "$error: min-qual must be an integer MAPQ\n";
    $minClip =~ m|\D| and die "$error: min-clip must be an integer number of bp\n";
    if ($bam) {
        if ($coordinate) {
            -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
            $tmpBam = "$tmpDir/svtools.$libKey.resort.".int(rand(1e6));
            #$stepsAfter = "samtools sort -no -m $maxMem - $tmpBam | ";
            $stepsAfter = "samtools sort -n -m $maxMem -T $tmpBam - | "; # new samtools sort format  
        }
        $stepsAfter .= "samtools view -h -";
        $repeat and $stepsBefore = "samtools view -hSb -";
    }
}

# main execution block
sub svtools_extract {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";

    # initialize
    setOptions_extract();    
    print STDERR "$utility extract: " . getTime(), "\n";
    initializeExclude($excludeFile);

    # open output streams
    my $tmpFile = getTmpFile("anomalies", $libKey, "gz");
    my $outFile = getOutFile("anomalies", $libKey, "gz");
    openOutputStream($tmpFile, \$tmpH, $TRUE);  
    openOutputStream($outFile, \$outH, $TRUE);
    $repeat and openOutputStream(undef, \$rptH, undef, undef, $stepsBefore);
    
    # do the work
    processInputFiles(\$inH, \@globs, $TRUE, $stepsAfter, $FALSE, \&handle_extract_file);
    closeHandles($rptH, $tmpH);    

    # commit unique read pairs from temporary anomalies file
    openInputStream($tmpFile, \$inH, $FALSE, "sort -S $maxMem"."b -k1,1 -k2,2");
    commitUniqAnom();
    closeHandles($inH, $outH);     
    unlink $tmpFile;

    # create histogram plots
    if($plotDir){
        print STDERR "$utility extract: creating size histograms\n";      
        createHist($tLens, 'tLens', $maxTLen, 'template length (bp)');
        createHist($rLens, 'rLens', $maxRLen, 'read length (nt)');
        createHist($anom{del}, 'del', $maxTLen*2, 'deletion size');
        createHist($anom{dup}, 'dup', $maxTLen*2, 'duplication size');
        createHist($anom{inv}, 'inv', $maxTLen*2, 'inversion size');   
    }

    # report counts
    reportCount($libType,    "$libKey\tlibrary type");
    reportCount($medRLen,    "$libKey\tmedian read length");
    reportCount($maxRLen,    "$libKey\tlargest read length");
    reportCount($medTLen,    "$libKey\tmedian expected insert size");
    reportCount($minTLen,    "$libKey\tsmallest expected insert size");
    reportCount($maxTLen,    "$libKey\tlargest expected insert size");
    reportCount($nPairs,     "$libKey\tread pairs processed");
    reportCount($nGaps,      "$libKey\tanomalous gaps between paired reads");
    reportCount($nSplits,    "$libKey\tsplit/chimeric reads");
    reportCount($nBadSplits, "$libKey\tbad reads with exactly two hard-clipped alignments");
    reportCount($nClips,     "$libKey\treads clips >= $minClip bases");
}

# process each input file
my $fileN = 1;
sub handle_extract_file {
    my ($file) = @_;
    
    # get pair statistics
    if ($fileN == 1) {
        print STDERR "$utility extract: getting read pair statistics\n";
        ($libType, $minTLen, $maxTLen, $medTLen, $maxRLen, $medRLen, $tLens, $rLens, $pTLen) =
            getStats($inH, $repeat, $rptH, \$statsLines, $nStatPairs, $minQual);  
        open(my $statsH, "<", \$statsLines) or die "$error: could not open handle to stats lines: $!\n";
        run_extract($statsH);
        close $statsH;
        $statsLines = ""; # release memory
        print $outH join("\t", "#application:svtools",
                                "version:$version",
                                "file_type:extract",
                                "sample:$sample",
                                "library:$library",
                                "libType:$libType",
                                "maxTLen:$maxTLen",
                                "medTLen:$medTLen",
                                "medRLen:$medRLen",
                                "minClip:$minClip",
                                "nStatPairs:$nStatPairs",
                                "pTLen:$pTLen"), "\n";     
    }

    # run anomalies extractions to tmp file
    print STDERR "$utility extract: $file\n";
    eof($inH) or run_extract($inH);
    closeHandles($inH);
    $tmpBam and -e $tmpBam and unlink $tmpBam;
    $fileN++;    
}

# extract anomalous read pairs
sub run_extract {
    my ($inH_) = @_;
    while (my $line = <$inH_>) {
        $line =~ m/^\@/ and next;
        $repeat and print $rptH $line;  # repeated output is NOT duplicate purged
        my ($qName) = split("\t", $line, 2);
        if($qName ne $prevName){
            processPair();
            $prevName = $qName;
            @alns = ();
        }
        push @alns, [split("\t", $line, 7)]; # last split = CIGAR (6, index 5)  
    }
    processPair(); # handle last pair
}
sub commitUniqAnom {
    my ($prevKey, $allowName) = ("NA");
    while (my $line = <$inH>) {
        my ($key, $name, $anom) = split("\t", $line, 3);
        if ($key ne $prevKey) {
            $prevKey = $key;
            $allowName = $name; # take all anomalies from _first_ read pair of a duplicate group
        }
        $name eq $allowName and print $outH $anom;        
    }
}

# parse paired alns that nominate a junction
sub processPair {
    @alns or return;
    $nPairs++;

    # enforce minimum mapping quality
    # ensure unique/non-duplicate DNA fragments by creating key based on all alns for read pair    
    # enforce reqion exclusions
    my ($maxAlnQ, @key) = (0);
    foreach my $aln(@alns){ # chrom, strand, pos, CIGAR defines an aligned segment
        $maxAlnQ >= $$aln[4] or $maxAlnQ = $$aln[4];
        my $pos; # the outermost position of the alignment IF all bases had mapped, including clips
        if($$aln[1] & 0x010){
            $pos = getEnd($$aln[3], $$aln[5]) + ($$aln[5] =~ m/$rgtClpRex/ ? $1 : 0);
        } else {
            $pos = $$aln[3] - ($$aln[5] =~ m/$lftClpRex/ ? $1 : 0);
        } # keep info on 1st/2nd read to allow reverse pairs to continue (they cannot be chimeric)
        push @key, join(":", $$aln[2], $$aln[1] & 0x010, $pos, $$aln[1] & 0x040);
        #push @key, join(":", $$aln[2], $$aln[1] & 0x010, $$aln[3], $$aln[5]);
        checkExclude($$aln[2], $$aln[3] - 1, $$aln[3]) and return; # just check mapped alignment position of read
    }                                                              # ignore pair if any read excluded  
    $maxAlnQ >= $minQual or return;  # at least one read had a quality > min-qual
    $key = join("::", sort @key); 

    # pair had no split reads
    if (@alns == 2) {
        if ($alns[0][1] & 0x002) { # a proper pair
            getOutInfo(\my @vs, undef, @alns);
            printOutClps(\@vs, @alns);
        } elsif($alns[0][1] & 0x00C){ # at least one unmapped read
            printAllClps(@alns);
        } else { # both mapped, not proper
            printGap($alns[0], undef, $alns[1], undef);
        }

    # at least one read had a split alignment
    } elsif(@alns > 2) {
        my (@alns1, @alns2);
        foreach my $aln(@alns){ # find the split read(s)
            $$aln[1] & 0x040 ? push @alns1, $aln : push @alns2, $aln;
        }
        (@alns1 > 2 or @alns2 > 2) and return; # for now, ignore highly complex reads
        printGap((@alns1 > 1 ? printSplit(\@alns1) : ($alns1[0], undef)),
                 (@alns2 > 1 ? printSplit(\@alns2) : ($alns2[0], undef)));
    }  
}

# collect basic information on reads, only as much as needed in steps
sub getOutInfo { # outermost pos and clip: outPos is the outermost MAPPED read position (excluding clips)
    my ($vs, $forcePos, @alns) = @_;
    my ($mRnm, $mPos, $tLen);
    foreach my $i(0, 1){
        my %vs;        
        $vs{str} = $alns[$i][1] & 0x010;
        if ($vs{str}) { # reverse
            $vs{outClp} = $alns[$i][5] =~ m/$rgtClpRex/ ? $1 : 0;
            ($forcePos or $vs{outClp} >= $minClip) and $vs{outPos} = getEnd($alns[$i][3], $alns[$i][5]);
            if ($vs{outClp} >= $minClip) {
                $vs{outSde} = "R";
                ($mRnm, $mPos, $tLen, $vs{seq}) = split("\t", $alns[$i][6], 5);
                $vs{outSeq} = substr($vs{seq}, -$vs{outClp});
            }
        } else { # forward
            $vs{outClp} = $alns[$i][5] =~ m/$lftClpRex/ ? $1 : 0;
            ($forcePos or $vs{outClp} >= $minClip) and $vs{outPos} = $alns[$i][3]; 
            if ($vs{outClp} >= $minClip) {
                $vs{outSde} = "L";
                ($mRnm, $mPos, $tLen, $vs{seq}) = split("\t", $alns[$i][6], 5);
                $vs{outSeq} = substr($vs{seq}, 0, $vs{outClp});
            }
        }
        $$vs[$i] = \%vs;
    }
}
sub getInnInfo { # innermost pos and clip
    my ($vs, @alns) = @_;
    my ($mRnm, $mPos, $tLen);
    foreach my $i(0, 1){
        my $vs = $$vs[$i];
        $$vs{seq} or ($mRnm, $mPos, $tLen, $$vs{seq}) = split("\t", $alns[$i][6], 5);
        if ($$vs{str}) { # reverse
            $$vs{innClp} = $alns[$i][5] =~ m/$lftClpRex/ ? $1 : 0 ;
            if ($$vs{innClp} >= $minClip) {
                $$vs{innPos} = $alns[$i][3]; 
                $$vs{innSde} = "L";
                $$vs{innSeq} = substr($$vs{seq}, 0, $$vs{innClp});
            }
        } else { # forward
            $$vs{innClp} = $alns[$i][5] =~ m/$rgtClpRex/ ? $1 : 0;  
            if ($$vs{innClp} >= $minClip) {
                $$vs{innPos} = getEnd($alns[$i][3], $alns[$i][5]);
                $$vs{innSde} = "R";
                $$vs{innSeq} = substr($$vs{seq}, -$$vs{innClp});
            }
        } 
    }
}

# handle the gap between reads
sub printGap {
    my ($aln1, $seq1, $aln2, $seq2) = @_;
    
    # suppress bad splits (see below)
    ($aln1 and $aln2) or return;

    # suppress unmapped partner of split read
    ($$aln1[1] & 0x004 or $$aln2[1] & 0x004) and return;

    # gather information
    getOutInfo(\my @vs, 1, $aln1, $aln2);

    # send insertions and expected to coverage map
    my $tLen = $vs[1]{outPos} - $vs[0]{outPos}; # can't trust tLen for split read
    if($$aln1[2] eq $$aln2[2] and
           ((!$vs[0]{str} and $vs[1]{str} and
              $tLen > 0 and  $tLen <= ($maxTLen - $vs[0]{outClp} - $vs[1]{outClp})) or
            (!$vs[1]{str} and $vs[0]{str} and
             -$tLen > 0 and -$tLen <= ($maxTLen - $vs[0]{outClp} - $vs[1]{outClp})))){
        # do nothing; handled by 'coverage' now
        
    # gather more information about anomalous pairs
    } else {
        $seq1 and $vs[0]{seq} = $seq1; # override hard clipped seqs
        $seq2 and $vs[1]{seq} = $seq2;          
        getInnInfo(\@vs, $aln1, $aln2);
        $vs[0]{maxDst} = $maxTLen - $vs[0]{outClp} - length($vs[1]{seq}) + $vs[1]{innClp};
        $vs[1]{maxDst} = $maxTLen - $vs[1]{outClp} - length($vs[0]{seq}) + $vs[0]{innClp};                      
        $vs[0]{minDst} = length($vs[0]{seq}) - $vs[0]{innClp} - $vs[0]{outClp};
        $vs[1]{minDst} = length($vs[1]{seq}) - $vs[1]{innClp} - $vs[1]{outClp};
        # strictly order the 1st and 2nd reads of an anomalous pair and print the gap
        if ($$aln1[2] eq $$aln2[2]) {
            my $size = abs($$aln2[3] - $$aln1[3]);
            if ($vs[0]{str} == $vs[1]{str}) { # order inversion by pos
                $size < $maxTLen*2 and $anom{inv}{$size}++;
                if ($$aln1[3] > $$aln2[3]) {
                    printGap_(\@vs, 1, 0, $aln2, $aln1);
                } else {
                    printGap_(\@vs, 0, 1, $aln1, $aln2)
                }
            } else {
                $size < $maxTLen*2 and $anom{
                    $vs[0]{str} ?
                    ($$aln1[3] > $$aln2[3] ? "del" : "dup"):
                    ($$aln2[3] > $$aln1[3] ? "del" : "dup")
                }{$size}++;
                if($vs[0]{str} > $vs[1]{str}) { # order del/dup by strand: 0,16
                    printGap_(\@vs, 1, 0, $aln2, $aln1);
                } else {
                    printGap_(\@vs, 0, 1, $aln1, $aln2)  
                }
            }
        } elsif($$aln1[2] gt $$aln2[2]){ # order trans by chrom
            printGap_(\@vs, 1, 0, $aln2, $aln1);
        } else {
            printGap_(\@vs, 0, 1, $aln1, $aln2)
        }        
    }

    # extract outer clips (suppress for split reads, already handled)
    $seq1 and $vs[0]{noOutClp} = 1;
    $seq2 and $vs[1]{noOutClp} = 1;
    printOutClps(\@vs, $aln1, $aln2);
}
sub printGap_ {
    my ($vs, $i1, $i2, $aln1, $aln2) = @_;
    $nGaps++;
    print $tmpH join("\t",  $key, $prevName, $_gap,
                            $$aln1[2], $$vs[$i1]{str}, $$aln2[2], $$vs[$i2]{str},
                            $$vs[$i1]{outPos}, $$vs[$i1]{maxDst}, $$vs[$i1]{minDst},
                            $$vs[$i2]{outPos}, $$vs[$i2]{maxDst}, $$vs[$i2]{minDst},
                            $$vs[$i1]{innSde} ?
                                ($$vs[$i1]{innPos}, $$vs[$i1]{innSde}, $$vs[$i1]{innSeq}, $$vs[$i1]{seq}) :
                                (0, "*", "*", "*"),
                            $$vs[$i2]{innSde} ?
                                ($$vs[$i2]{innPos}, $$vs[$i2]{innSde}, $$vs[$i2]{innSeq}, $$vs[$i2]{seq}) :
                                (0, "*", "*", "*")), "\n";  
}

# handle split reads, i.e. sequenced junctions
sub printSplit {
    my ($alns) = @_;

    # get the full read sequence (i.e not hard clipped)
    my $k = $$alns[0][5] =~ m/H/ ? 1 : 0;
    if($k == 1 and $$alns[$k][5] =~ m/H/){
        $nBadSplits++;   # two hard-clipped reads, can't get the sequence!
        return (0, 0);   # not clear why this _rarely_ happens since it is filtered to two max alns 
    }                      
    my ($mRnm, $mPos, $tLen, $seq) = split("\t", $$alns[$k][6], 5);
    if ($k == 1 and ($$alns[0][1] & 0x010) != ($$alns[1][1] & 0x010)) {
        $seq =~ tr/ACGT/TGCA/;
        $seq = reverse($seq);  
    }
    my $readLen = length($seq);

    # characterize the halves of the junction
    my ($i1, $i2, @vs) = (0, 1);    
    foreach my $i($i1, $i2){
        my %vs;
        $vs{end} = getEnd($$alns[$i][3], $$alns[$i][5]);    
        $vs{lftClp} = $$alns[$i][5] =~ m/$lftClpRex/ ? $1 : 0;
        $vs{rgtClp} = $$alns[$i][5] =~ m/$rgtClpRex/ ? $1 : 0;
        ($vs{jxnSde}, $vs{jxnClp}, $vs{outClp}, $vs{jxnPos}) =
            $vs{lftClp} > $vs{rgtClp} ?
                ("L", $vs{lftClp}, $vs{rgtClp}, $$alns[$i][3]) :
                ("R", $vs{rgtClp}, $vs{lftClp}, $vs{end});
        $vs{maxDst} = $readLen - $vs{outClp} - $vs{jxnClp}; 
        $vs{str} = $$alns[$i][1] & 0x010;
        $vs{invStr} = (( $vs{str} and $vs{jxnSde} eq "R") or 
                       (!$vs{str} and $vs{jxnSde} eq "L")) ?
                            ($vs{str} ? 0 : 16):
                             $vs{str}; # for pairing with gaps, invert the innermost alignment
        $vs{idxPos} = $vs{invStr} ? $vs{end} : $$alns[$i][3];
        $vs[$i] = \%vs;
    }

    # strictly order the halves
    if ($$alns[$i1][2] eq $$alns[$i2][2]) {
        if ($vs[$i1]{invStr} == $vs[$i2]{invStr}) { # order inversion by pos
            if ($$alns[$i1][3] > $$alns[$i2][3]) {
                ($i1, $i2) = invertJxn($i1, $i2, \@vs, \$seq); 
            }
        } elsif($vs[$i1]{invStr} > $vs[$i2]{invStr}) { # order del/dup by strand
            ($i1, $i2) = invertJxn($i1, $i2, \@vs, \$seq); 
        }
    } elsif($$alns[$i1][2] gt $$alns[$i2][2]){ # order trans by chrom
        ($i1, $i2) = invertJxn($i1, $i2, \@vs, \$seq); 
    }
    sub invertJxn {
        my ($i1, $i2, $vs, $seq) = @_;
        ($i1, $i2) = ($i2, $i1);
        unless($$vs[$i1]{str} == $$vs[$i2]{str}){
            $$seq =~ tr/ACGT/TGCA/;
            $$seq = reverse($$seq);    
        }
        return ($i1, $i2);
    }

    # characterize the junction
    my $microH = $readLen + $vs[$i1]{outClp} + $vs[$i2]{outClp} -
                 $vs[$i1]{lftClp} - $vs[$i1]{rgtClp} - $vs[$i2]{lftClp} - $vs[$i2]{rgtClp}; 
    my $insSeq = $microH < 0 ?  # negative microH = length of insSeq
                    ($vs[$i1]{jxnSde} eq "R" ?
                        substr(substr($seq, -$vs[$i1]{rgtClp}), 0, -$microH) :
                        substr(substr($seq, 0, $vs[$i1]{lftClp}), $microH) ) :
                    "-";
    $microH < 0 and $microH = 0;

    # print the split junction
    $nSplits++;
    print $tmpH join("\t",  $key, $prevName, $_split,
                            $$alns[$i1][2], $vs[$i1]{invStr}, $$alns[$i2][2], $vs[$i2]{invStr},
                            $vs[$i1]{idxPos}, $vs[$i1]{maxDst}, $vs[$i1]{maxDst},
                            $vs[$i2]{idxPos}, $vs[$i2]{maxDst}, $vs[$i2]{maxDst},
                            $vs[$i1]{jxnPos}, $vs[$i1]{jxnSde},
                            $vs[$i2]{jxnPos}, $vs[$i2]{jxnSde},
                            $microH, $insSeq, $seq), "\n";

    # return innermost alignment for calling the gap
    my $j = $i1;
    unless((  $vs[$j]{str} and $vs[$j]{jxnSde} eq "R") or
            (!$vs[$j]{str} and $vs[$j]{jxnSde} eq "L")){
        $j = $i2;
        unless($vs[$i1]{str} == $vs[$i2]{str}){
            $seq =~ tr/ACGT/TGCA/;
            $seq = reverse($seq);    
        } 
    }
    return ($$alns[$j], $seq); 
}

# print large clips as ancillary junction support
sub printAllClps { # for reads with an unmapped partner
    foreach my $aln(@_){
        unless($$aln[1] & 0x004){
            if ($$aln[5] =~ m/$lftClpRex/ and $1 >= $minClip) {
                my ($mRnm, $mPos, $tLen, $seq) = split("\t", $$aln[6], 5);
                $nClips++;
                print $tmpH join("\t",  $key, $prevName, $_clip,
                                        $$aln[2], $$aln[3], "L",
                                        substr($seq, 0, $1), $seq), "\n";    
            }
            if ($$aln[5] =~ m/$rgtClpRex/ and $1 >= $minClip) {
                my ($mRnm, $mPos, $tLen, $seq) = split("\t", $$aln[6], 5);
                $nClips++;
                print $tmpH join("\t",  $key, $prevName, $_clip,
                                        $$aln[2], getEnd($$aln[3], $$aln[5]), "R",
                                        substr($seq, -$1), $seq), "\n";    
            }
        }
    }  
}
sub printOutClps { # the outer clips of a mapped pair (inner clips handled above)
    my ($vs, @alns) = @_;
    foreach my $i(0, 1){
        (!$$vs[$i]{outSde} or $$vs[$i]{noOutClp}) and next;
        $nClips++;
        print $tmpH join("\t",  $key, $prevName, $_clip,
                                $alns[$i][2], $$vs[$i]{outPos}, $$vs[$i]{outSde},
                                $$vs[$i]{outSeq}, $$vs[$i]{seq}), "\n"; 
    }
}


# create histogram plots
sub createHist {
    my ($hash, $fileType, $maxX, $xLab) = @_;
    my $imgFile = "$plotDir/$libKey.$fileType.jpg";
    open(my $outH, "|-", "Rscript $libDir/histogram.R $imgFile $libKey $maxX '$xLab'")
        or die "$error: createHist could not open output stream for writing: $!\n"; 
    foreach my $key(sort {$a <=> $b} keys %$hash){
        print $outH join("\t", $key, $$hash{$key}), "\n";
    }
    close $outH;
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
#7	MRNM	Mate Reference sequence NaMe (= if same as RNAME)
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
