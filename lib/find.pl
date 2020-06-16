use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            $_gap $_split $_clip $pTLenBin
            %idxSubs $idxH $libKey);
require "$libDir/index_common.pl";
my (@globs, $group, $genomeFasta, $maxMem, 
    $minTLen, $minSetPairs, $unique, $uniqueStrict, $jxnOnly, $minSize, $maxSize,
    $chrom, $chroms, $useMask);
my ($gZip, $stepsAfter, $stepsBefore);
my $minMinClip = 1e9;
my $maxMaxTLen = 0;
my $wrkMaxTLen = 0;
my (%chroms, @chroms, %samples, @samples, %libKeys, %extDirs, %clips, %mask);
my ($chr1, $str1, $chr2, $str2);
my (@pairs, @ks, %svSigs, %evid);
my $svID = 0;
my ($nClips, $nPairs, $nPairsUsed, $nFindAttempts,
    $nReconstrs, $nSplits, $nGaps,
    $nDupSets, $nSets, $nUnique, $nUniqueStrict) = (0) x 100;

# manage options
sub setOptions_find {
    setOptionValue(\$group,         'group');
    setOptionValue(\$genomeFasta,   'genome-fasta');
    setOptionValue(\$outputDir,     'output-dir');  
    setOptionValue(\$tmpDir,        'tmp-dir',       '/tmp');
    setOptionValue(\$maxMem,        'max-mem',       1000000000);
    setOptionValue(\$minTLen,       'min-tLen',      '1.5/-1/0');
    setOptionValue(\$minSetPairs,   'min-set-pairs', 4);
    setOptionValue(\$unique,        'unique');
    setOptionValue(\$uniqueStrict,  'unique-strict');
    setOptionValue(\$jxnOnly,       'jxn-only');
    setOptionValue(\$minSize,       'min-size',      0);
    setOptionValue(\$maxSize,       'max-size',      1000000000);
    setOptionValue(\$chrom,         'chrom');
    setOptionValue(\$chroms,        'chroms');
    setOptionValue(\$useMask,       'use-mask');
    #-------------------------------------
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    -e $genomeFasta or die "$error: file not found: $genomeFasta\n";
    my $genIdx = "$genomeFasta.fai";
    -e $genIdx or system("samtools faidx $genomeFasta");
    $maxMem =~ m|\D| and die "$error: max-mem must be an integer number of bytes\n";
    $minSetPairs =~ m|\D| and die "$error: min-set-pairs must be an integer number\n";
    $minSize =~ m|\D| and die "$error: min-size must be an integer number\n";
    $maxSize =~ m|\D| and die "$error: max-size must be an integer number\n";
    %chroms = map { $_ => 1} ($chrom, split(",", $chroms));
    @chroms = sort {$a cmp $b} keys %chroms;
}

# main execution block
sub svtools_find {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";
    
    # initialize
    setOptions_find();    
    print STDERR "$utility find: " . getTime(), "\n";
    
    # pull the appropriate data from the set of input extract files
    print STDERR "$utility find: assembling $chrom anomalous pair information\n";
    my $tmpFile = getTmpFile("pairs", "$group.$chrom", "gz");
    openOutputStream($tmpFile, \$tmpH, $TRUE);
    processInputFiles(\$inH, \@globs, $TRUE, $FALSE, $TRUE, \&handle_find_file);
    closeHandles($tmpH);
    @samples = sort {$a cmp $b} keys %samples;

    # find sets in previously extracted data, by chrom/strand combination
    print STDERR "$utility find: finding sets\n";
    openInputStream($tmpFile, \$inH, $FALSE, "sort -S $maxMem"."b -k1,4");
    #my $smpFile = getOutFile("samples", "$group.$chrom", "txt");
    #my $outFile = getOutFile("sets",    "$group.$chrom", "bgz");
    #openOutputStream($smpFile, \my $smpH);
    #openOutputStream($outFile, \$outH, undef, undef,
    #                 "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    my $fileID = "$group.$chrom";
    my $outFile = getOutFile("sets", $fileID, "gz");
    openOutputStream($outFile, \$outH, $TRUE);
    my $idxFile = getOutFile('anomalies', $fileID, "gz");
    openOutputStream($idxFile, \$idxH, $TRUE);
    print $outH join("\t", "#application:svtools",
                            "version:$version",
                            "file_type:find",
                            "group:$group",
                            "minTLen:$minTLen",
                            "minSetPairs:$minSetPairs",
                            "minSize:$minSize",
                            "maxSize:$maxSize",
                            "chroms:".join(",", @chroms),
                            "samples:".join(",", @samples)), "\n";      
    threadInput();
    #print $smpH join("\t", @samples), "\n";
    closeHandles($inH, $outH, $idxH);    
    unlink $tmpFile;
    #system("tabix -p bed $outFile");

    # report counts
    print STDERR "$utility find: $group: ", join(" ", sort {$a cmp $b} keys %libKeys), "\n";
    reportCount($maxMaxTLen,    "largest expected insert size over all samples");
    reportCount($wrkMaxTLen,    "bp flanking reported junctions");
    reportCount($minMinClip,    "smallest allowed clip length over all samples");
    reportCount($nClips,        "outer clipped reads across all samples");
    reportCount($nPairs,        "anomalous read pairs across all samples");
    reportCount($nPairsUsed,    "anomalous read pairs >min-tLen and used in find");
    reportCount($nFindAttempts, "attempts to find a structural variant junction");
    reportCount($nReconstrs,    "candidate junctions reconstructed from inner clips");
    reportCount($nSplits,       "candidate junctions obtained from split reads");
    reportCount($nGaps,         "candidate structural variants with only read pair evidence");
    reportCount($nDupSets,      "structural variants disregarded as duplicates");    
    reportCount($nSets,         "called structural variants");
    reportCount($nUnique,       "structural variants called in only one sample");
    reportCount($nUniqueStrict, "structural variants present in only one sample");
}

# load clips while discovering sample information
sub handle_find_file {
    my ($file) = @_;
    
    # read the extract header
    my (%minTLen);
    my $line = <$inH>;
    chomp $line;
    $line =~ m/^#(.+)/ or die "$utility find error: no header line in extract file: $file\n";       
    my %tags = getFileHeader($1, 'extract', qw(sample library pTLen nStatPairs
                                               medRLen medTLen maxTLen minClip));
    my $libKey  = "$tags{sample}.$tags{library}";
    $libKeys{$libKey} and die "$error: library $libKey encountered >1 time in inputs\n";
    $file =~ m/(.+)\//g;
    $extDirs{$libKey} =  $1;
    my %pTLen = map {
        my ($key, $val) = $_ =~ m/(.+)=(.+)/ ? ($1, $2) : ('xxx', 0);
        $key => $val;
    } split(";", $tags{pTLen});
    $pTLen{extreme} = 1/$tags{nStatPairs};            
    $libKeys{$libKey} = [$tags{medRLen}, $tags{medTLen}, $tags{maxTLen}, $tags{minClip}, \%pTLen];
    setMinTLen('find', $libKey, $minTLen, \%minTLen, \%tags);
    push @{$samples{$tags{sample}}}, $tags{library};
    $maxMaxTLen >= $tags{maxTLen} or $maxMaxTLen = $tags{maxTLen};
    $minMinClip <= $tags{minClip} or $minMinClip = $tags{minClip};   

    # pull the chromosome-specific pairs
    while (my $line = <$inH>) {
        chomp $line;
        # restrict attention to the index chromosome (including all of its sorted partners (even some we may not want))
        my ($pairType, $rest) = split("\t", $line, 2);
        $rest =~ m/^$chrom\t/ or next;
        # save the (mostly outer) clips in memory  NOTE: takes a lot of memory!!
        if ($pairType == $_clip) {
            $nClips++;
            #my ($chr, $jxnPos, $jxnSde, $clpSeq, $seq) = split("\t", $rest); 
            #push @{$clips{join(",", $chr, $jxnPos, $jxnSde)}}, [$clpSeq, $seq, $libKey];
        # save gaps and splits on disk
        } else {
            $nPairs++;
            checkTLen($pairType, $rest, \%minTLen) or next;
            $nPairsUsed++;
            print $tmpH join("\t", $rest, $tags{sample}, $libKey, $pairType), "\n";
        }               
    }
    closeHandles($inH);
    $wrkMaxTLen = int($maxMaxTLen/100+0.99999)*100; # round jxn padding up to 100  
}

# process data by pairs of chr1:str1 with chr2:str2
sub threadInput {
    @pairs = ();
    my $prevKey;
    while (my $line = <$inH>) {
        chomp $line;
        my ($chr1_, $str1_, $chr2_, $str2_, @rest) = split("\t", $line);
        $chroms{$chr2_} or next; # skip partners we don't want
        if (!$prevKey){
            $prevKey = "$chr1_,$str1_,$chr2_,$str2_";
        } elsif($prevKey ne "$chr1_,$str1_,$chr2_,$str2_"){
            ($chr1, $str1, $chr2, $str2) = split(",", $prevKey);
            print STDERR "$utility find: ".join("\t", $chr1, $str1, $chr2, $str2)."\n";
            @pairs >= $minSetPairs and threadP1();
            @pairs = ();
            $prevKey = "$chr1_,$str1_,$chr2_,$str2_";
        }
        push @pairs, \@rest;  
    }
    ($chr1, $str1, $chr2, $str2) = split(",", $prevKey);
    @pairs >= $minSetPairs and threadP1();
}

# find definitive hard breaks in P1 against maxMaxTlen criterion by linear threading
sub threadP1 {
    loadMask($chr1, $chr2);
    my @parseParam =
        (0, 1, 2,
        \&parsePMatrix, 3, 4, 5,
        \&processNoJunction);
    my ($prevPair, @is);
    foreach my $i(sort {$pairs[$a][0] <=> $pairs[$b][0]} 0..$#pairs){
        checkMask($pairs[$i]) or next;
        if ($prevPair and ($pairs[$i][0] - $$prevPair[0]) > $maxMaxTLen){
            @is >= $minSetPairs and parsePMatrix(\@is, @parseParam);
            @is = ();  
        }
        push @is, $i;
        $prevPair = $pairs[$i];
    }
    @is >= $minSetPairs and parsePMatrix(\@is, @parseParam);
}

# load the mask list of read alignments that exist on the ends of observed proper pairs
sub loadMask {
    $useMask or return;
    return;
    # masking is awaiting the availability of UMIs
    #my (@chroms) = @_;
    #my %chroms = map { $_ => 1 } @chroms;
    #foreach my $chrom(keys %mask){ # delete unused masks for memory efficiency
    #    $chroms{$chrom} or delete $mask{$chrom};
    #}
    #foreach my $chrom(@chroms){
    #    $mask{$chrom} and next; # continuing with a previous chromosome
    #    foreach my $libKey(keys %libKeys){
    #        $mask{$chrom}{$libKey} = {};
    #        $inputDir = $extDirs{$libKey};
    #        my $maskFile = getInFile('mask', "knowns", $libKey, "bgz");
    #        my $stream = "tabix $maskFile $chrom";
    #        open my $mskH, "-|", $stream or die "$error: could not open $maskFile\n";
    #        while(my $line = <$mskH>){
    #            chomp $line;
    #            my ($chr, $pos, $str) = split("\t", $line); # key = out-pos, str(+/-) [, readN(1/2)]
    #            $mask{$chrom}{$libKey}{join(":", $pos, $str)}++;
    #        }
    #        closeHandles($mskH);
    #    }
    #}  
}
sub checkMask {
    $useMask or return 1;
    return 1;
    # masking is awaiting the availability of UMIs
    #my ($pair) = @_;
    #$mask{$chr1}{$$pair[$#$pair-1]}{join(":", $$pair[0], $str1)} and return;
    #$mask{$chr2}{$$pair[$#$pair-1]}{join(":", $$pair[3], $str2)} and return;
    #return 1;
}

# solve for all largest unique non-breaking pos sets from current sub-matrix of pairs
sub parsePMatrix {
    my ($is, $pos, $maxDst, $minDst, $nextSub, @nextParam) = @_;
    my (@sets, %seen, %sent);
    foreach my $i1(@$is){
        my @set;
        foreach my $i2(@$is){
            abs($pairs[$i2][$pos] - $pairs[$i1][$pos]) <=
                ($str1 ?
                 $pairs[$i2][$maxDst] - $pairs[$i1][$minDst] :
                 $pairs[$i1][$maxDst] - $pairs[$i2][$minDst])
            and push @set, $i2;
        }
        @set >= $minSetPairs and push @sets, \@set;
    }   
    P_SET: foreach my $set(sort {@$b <=> @$a} @sets){ # largest sets first
        my $id = join(",", sort {$a <=> $b} @$set);
        $seen{$id} and next;
        $seen{$id}++;
        foreach my $sent(keys %sent){
            my (@common, %n);
            foreach my $i(@$set, split(",", $sent)) { $n{$i}++ }
            foreach my $i(keys %n){ $n{$i} > 1 and push @common, $i }
            $id eq join(",", sort {$a <=> $b} @common) and next P_SET; # subset of sent/larger set
        }
        &$nextSub($set, @nextParam);
        $sent{$id}++;
    }
}

# determine whether split and clipped reads establish a definitive junction
sub findJunction {
    $nFindAttempts++;
     
    #######################
    processP2Group();
    return;

    # collect evidence from split reads and inner clips (present in pairs set)
    my (%splits, %spltInf, %pairClps, %clpInf);    
    foreach my $k(@ks){        
        if ($pairs[$k][$#{$pairs[$k]}] == $_split) {
            my $clp1 = join(",", $chr1, @{$pairs[$k]}[6..7]);
            my $clp2 = join(",", $chr2, @{$pairs[$k]}[8..9]);
            $splits{"$clp1,$clp2"}++;
            $spltInf{"$clp1,$clp2"}{join(",", @{$pairs[$k]}[10,11])}++;
            $pairClps{$clp1}++;
            $pairClps{$clp2}++; # do not have or need clp info, it is in split
        }
        #else {
        #    if ($pairs[$k][6]) {
        #        my $clp = join(",", $chr1, @{$pairs[$k]}[6..7]);
        #        $pairClps{$clp}++;
        #        push @{$clpInf{$clp}}, [@{$pairs[$k]}[8..9,$#{$pairs[$k]}-1], $k];
        #    }
        #    if ($pairs[$k][10]) {
        #        my $clp = join(",", $chr2, @{$pairs[$k]}[10..11]);
        #        $pairClps{$clp}++;
        #        push @{$clpInf{$clp}}, [@{$pairs[$k]}[12..13,$#{$pairs[$k]}-1], $k];
        #    }
        #}
    }
    
    ## add evidence from matching outer clips (not present in pairs set)
    #foreach my $clp(keys %pairClps){
    #    if($clips{$clp}){
    #        $pairClps{$clp} += @{$clips{$clp}};
    #        push @{$clpInf{$clp}}, @{$clips{$clp}};
    #    }
    #}

    # find the most abundant clips, especially if they matched a split
    if (keys %pairClps >= 2) {
        my @pairClps = sort {$pairClps{$b} <=> $pairClps{$a}} keys %pairClps;
        
        # best clips match a split, obvious high quality junction, use it
        if ($splits{join(",", @pairClps[0..1])}){
            processJunction(join(",", @pairClps[0..1]), \%spltInf);
        } elsif($splits{join(",", reverse @pairClps[0..1])}) {
            processJunction(join(",", reverse @pairClps[0..1]), \%spltInf); 
            
        # clips only, still pretty good evidence
        } elsif(keys %splits == 0) {
            #searchClips(\%clpInf, \@pairClps);
            $jxnOnly or processP2Group();

        # clips don't match existing splits, work from best split instead  
        } else { 
            processJunction((sort {$splits{$b} <=> $splits{$a}} keys %splits)[0], \%spltInf);
        }
        
    # only one clip, try to turn it into a junction
    #} elsif (keys %pairClps == 1) {
        #searchClips(\%clpInf, [keys %pairClps]);
        
    # no clips to define a junction, proceed using VAMP logic  
    } else { 
        $jxnOnly or processP2Group();
    }
}

# try to find an inner clip on the other side of that pair to make a junction
sub searchClips {
    my ($clpInf, $pairClps) = @_;

    # collect information on the candidate inner clip
    foreach my $clp(@$pairClps){
        my ($jxnChr, $jxnPos, $jxnSde) = split(",", $clp);
        my %othClps;

        # sort the pairs containing this clip by clip length
        foreach my $clpI(sort {length($$b[0]) <=> length($$a[0])} @{$$clpInf{$clp}}){
            $$clpI[3] or next;  # k is undef if outer clip
            my ($othChr, $othJxnPos, $othJxnSde, $microH, $insSeq);                 

            # find the read with the inner clip and search the opposite side
            if ($jxnChr eq $chr1 and       # read 1 has the inner clip
                $jxnPos == $pairs[$$clpI[3]][6] and
                $jxnSde eq $pairs[$$clpI[3]][7]) {
                ($othChr, $othJxnPos, $othJxnSde, $microH, $insSeq) =
                    searchClip($clpI, $jxnSde, $chr2, $str2, 3..5);
            } elsif ($jxnChr eq $chr2 and  # read 1 has the inner clip
                     $jxnPos == $pairs[$$clpI[3]][10] and
                     $jxnSde eq $pairs[$$clpI[3]][11]) {
                ($othChr, $othJxnPos, $othJxnSde, $microH, $insSeq) =
                    searchClip($clpI, $jxnSde, $chr1, $str1, 0..2);
            }
            
            # other clip already known, high quality junction, use it
            if ($othChr){
                if($$clpInf{join(",", $othChr, $othJxnPos, $othJxnSde)}) {
                    return reconstructJunction($microH, $insSeq,
                                               [$jxnChr, $jxnPos, $jxnSde],
                                               [$othChr, $othJxnPos, $othJxnSde]);
                } else {
                    $othClps{join(",", $othChr, $othJxnPos, $othJxnSde, $microH, $insSeq)}++;
                }
            }
        }
        
        # choose the best other clip observed over all pairs
        if (keys %othClps) {
            my ($othChr, $othJxnPos, $othJxnSde, $microH, $insSeq) =
                split(",", (sort {$othClps{$b} <=> $othClps{$a}} keys %othClps)[0]);
                return reconstructJunction($microH, $insSeq,
                                           [$jxnChr, $jxnPos, $jxnSde],
                                           [$othChr, $othJxnPos, $othJxnSde]);
        }
    }
    
    # no other clip found, i.e. no junction to use
    $jxnOnly or processP2Group();            
}
sub searchClip { 
    my ($clpInf, $jxnSde, $othChr, $othStr, $outPosI, $maxDstI, $minDstI) = @_;
    my $k = $$clpInf[3];  # clpInf = [clpSeq, seq, sample, k]
    my $seq = $str1 == $str2 ? revComp($$clpInf[1]) : $$clpInf[1];   
    my $clpSeqLen = length($$clpInf[0]);

    # collect target sequence on other side of pair
    my $othSeq = getRegSeqs($othStr ?
                   "$othChr:".($pairs[$k][$outPosI]-$pairs[$k][$maxDstI])."-".
                              ($pairs[$k][$outPosI]-$pairs[$k][$minDstI]+$clpSeqLen) :
                   "$othChr:".($pairs[$k][$outPosI]+$pairs[$k][$minDstI]-$clpSeqLen)."-".
                              ($pairs[$k][$outPosI]+$pairs[$k][$maxDstI]));
    $othSeq = $$othSeq[0];
    my ($othJxnPos, $othJxnSde, $aln, $insSeq, $wrkClp, $microH);
    
    # try various insertion lengths to match other side
    foreach my $insSeqL(0..$clpSeqLen-$samples{$$clpInf[2]}[3]){
        my $wrkClpL = $clpSeqLen-$insSeqL;
        if (($str1 == $str2 and $jxnSde eq "R") or
            ($str1 != $str2 and $jxnSde eq "L")) {
            ($wrkClp, $insSeq, $aln) = $seq =~ m|(.{$wrkClpL})(.{$insSeqL})(.+)|;
        } else {
            ($aln, $insSeq, $wrkClp) = $seq =~ m|(.+)(.{$insSeqL})(.{$wrkClpL})|;
        }

        # clpSeq is found exactly one time on other side of junction?
        my $othPos = index($othSeq, $wrkClp);
        if ($othPos >= 0 and index(substr($othSeq, $othPos + 1), $wrkClp) == -1) {
            
            # if no insertion, look for a microhomology
            if($insSeqL){
                $microH = 0;
            } elsif ($othStr) {
                $microH = 1;
                while ($othPos and substr($aln,    -$microH,  1) eq
                                   substr($othSeq, $othPos-1, 1)) {
                    $microH++;
                    $othPos--;
                }
                $microH--;
            } else {
                $microH = 0;
                while (substr($aln,    $microH,  1) eq
                       substr($othSeq, $othPos+$clpSeqLen, 1)) {
                    $microH++;
                    $othPos++;
                } 
            }

            # name the junction
            if ($othStr) {
                $othJxnPos = $pairs[$k][$outPosI]-$pairs[$k][$maxDstI] + $othPos;
                $othJxnSde = "L";
            } else {
                $othJxnPos = $pairs[$k][$outPosI]+$pairs[$k][$minDstI] + $othPos - 1;
                $othJxnSde = "R";                    
            }
        }
        
        # this pair matched >= 1 time, shorter clips won't help if was promiscuous match
        $othPos >= 0 and last;
    }
    
    # return a found junction
    if ($insSeq and $str1 == $str2) { $insSeq = revComp($insSeq) }
    $othJxnSde and return ($othChr, $othJxnPos, $othJxnSde, $microH, $insSeq);
}

# build putative junction from _two_ clips when no split read gave it to us
sub reconstructJunction {
    my ($microH, $insSeq, @clps) = @_;

    # strictly order the clips
    my ($i1, $i2) = (0, 1);
    if ($clps[0][0] eq $clps[1][0]) {
        if ($clps[0][2] eq $clps[1][2]) { # order inversion by pos
            if ($clps[0][1] > $clps[1][1]) {
                ($i1, $i2) = (1, 0);
                $insSeq = revComp($insSeq);
            }
        } elsif($clps[0][2] lt $clps[1][2]) { # order del/dup by strand (L/R sort opposite 0/16)
            ($i1, $i2) = (1, 0);
        }
    } elsif($clps[0][0] gt $clps[1][0]){ # order trans by chrom
        ($i1, $i2) = (1, 0);
        $clps[0][2] eq $clps[1][2] and $insSeq = revComp($insSeq);
    } 

    # continue with the reassembled junction
    $nReconstrs++;
    my $jxn = join(",", @{$clps[$i1]}, @{$clps[$i2]});
    $insSeq or $insSeq = "-";
    processJunction($jxn, {$jxn => {join(",", $microH, $insSeq) => 1}}, 1);
}

# long reads have completely characterized an SV junction
sub processJunction {
    my ($jxn, $spltInf, $isClips) = @_;
    my ($chr1, $jxnPos1, $jxnSde1, $chr2, $jxnPos2, $jxnSde2) = split(",", $jxn);
    my $clp1 = join(",", $chr1, $jxnPos1, $jxnSde1);
    my $clp2 = join(",", $chr2, $jxnPos2, $jxnSde2);
    my ($microH, $insSeq) = # the most common jxn info
        split(",", (sort {$$spltInf{$jxn}{$b} <=> $$spltInf{$jxn}{$a}} keys %{$$spltInf{$jxn}})[0]);

    # assemble the junction sequence
    my $jxnSource = $isClips ? "clip" : "split";
    $isClips ? $nClips++ : $nSplits++;
    my $jxnSeq;
    
    #### MOST likeely a range error here (i.e. <0)
    my $refRegs = [$jxnSde1 eq "L" ?
                "$chr1:".$jxnPos1."-".($jxnPos1+$wrkMaxTLen) :
                "$chr1:".($jxnPos1-$wrkMaxTLen)."-".$jxnPos1,
                $jxnSde2 eq "L" ?
                "$chr2:".$jxnPos2."-".($jxnPos2+$wrkMaxTLen) :
                "$chr2:".($jxnPos2-$wrkMaxTLen)."-".$jxnPos2];
    my $flnkSeqs = getRegSeqs(@$refRegs);
    if ($jxnSde1 eq "R" and $jxnSde2 eq "L") {
        $jxnSeq = $$flnkSeqs[0].($insSeq eq "-" ? "" : $insSeq).substr($$flnkSeqs[1], $microH)
    } elsif($jxnSde1 eq "L" and $jxnSde2 eq "R") {
        $jxnSeq = $$flnkSeqs[1].($insSeq eq "-" ? "" : $insSeq).substr($$flnkSeqs[0], $microH)
    } elsif($jxnSde1 eq "R" and $jxnSde2 eq "R"){
        $jxnSeq = $$flnkSeqs[0].($insSeq eq "-" ? "" : $insSeq).substr(revComp($$flnkSeqs[1]), $microH)        
    } else {
        $jxnSeq = revComp(($insSeq eq "-" ? "" : $insSeq).substr($$flnkSeqs[0], $microH)).$$flnkSeqs[1]  
    }  

    # determine the junction type and event size
    my ($svSize, $jxnType);
    if ($chr1 ne $chr2) {
        $svSize = 0;        
        $jxnType = "trans";
    } elsif($jxnSde1 eq $jxnSde2){
        $svSize = $jxnPos2 - $jxnPos1 - $microH;         
        $jxnType = $jxnSde1 eq "R" ? "invF" : "invR";
    } else {
        $svSize = $jxnPos1 - $jxnPos2 + 1 - $microH;
        $jxnType = $svSize < 0 ? "del" : "dup";  # del svSize is a negative number
    }
    
    #########################
    return;
    

    # separate the read pairs that (do not) match this junction
    my ($maxEvid, @jxnKs, @nonJxnKs) = (0);
    initEvid();
    foreach my $k(@ks){
        #print join("\t", @{$pairs[$k]}[0,3]), "\n";
        # splits
        my $evidI;
        if ($pairs[$k][$#{$pairs[$k]}] == $_split) {
            $jxn eq join(",", $chr1, @{$pairs[$k]}[6..7], $chr2, @{$pairs[$k]}[8..9]) and $evidI = 1;
            
        # does read pair bracket the junction?
        } elsif (($jxnSde1 eq "R" ? $pairs[$k][0] <= $jxnPos1 : $pairs[$k][0] >= $jxnPos1) and
                 ($jxnSde2 eq "R" ? $pairs[$k][3] <= $jxnPos2 : $pairs[$k][3] >= $jxnPos2)){
            
            # if a read crosses the junction, is it inner-clipped appropriately?
            if ((length($pairs[$k][9]) - abs($pairs[$k][0] - $jxnPos1) - 1) >=
                    $libKeys{$pairs[$k][$#{$pairs[$k]}-1]}[3]) {
                $pairs[$k][6] and $clp1 eq join(",", $chr1, @{$pairs[$k]}[6..7]) and $evidI = 3;
            } elsif ((length($pairs[$k][13]) - abs($pairs[$k][3] - $jxnPos2) - 1) >=
                    $libKeys{$pairs[$k][$#{$pairs[$k]}-1]}[3]) {
                $pairs[$k][10] and $clp2 eq join(",", $chr2, @{$pairs[$k]}[10..11]) and $evidI = 5;
                
            # does junction predict an appropriate tLen for a gap fragment?
            } elsif(abs($jxnPos1 - $pairs[$k][0]) +
                    abs($jxnPos2 - $pairs[$k][3])
                        <= $libKeys{$pairs[$k][$#{$pairs[$k]}-1]}[2]){
                    $evidI = 4;
            }
        }
        if (defined $evidI) {
            push @jxnKs, $k;
            $evid{all}[0]++;
            $evid{all}[$evidI]++;
            $evid{$pairs[$k][$#{$pairs[$k]}-2]}[0]++;
            $evid{$pairs[$k][$#{$pairs[$k]}-2]}[$evidI]++;
            $maxEvid >= $evid{$pairs[$k][$#{$pairs[$k]}-2]}[0] or # best per-sample evidence
                $maxEvid = $evid{$pairs[$k][$#{$pairs[$k]}-2]}[0];
        } else {
            push @nonJxnKs, $k;
        }
    }

    # if this juncion has sufficient evidence for at least one sample, commit it
    if($maxEvid >= $minSetPairs){
        # add outer clip evidence (not counted against threshold since one-sided)
        if ($clips{$clp1}) {
            $evid{all}[2] = @{$clips{$clp1}};
            foreach my $clp(@{$clips{$clp1}}){ $evid{$$clp[2]}[2]++ }
        }
        if ($clips{$clp2}) {
            $evid{all}[6] = @{$clips{$clp2}};
            foreach my $clp(@{$clips{$clp2}}){ $evid{$$clp[2]}[6]++ }
        } 
        commitSV(\@jxnKs, $jxnSource, $jxnType, $svSize, 0, 0, 0, 0, # TODO: needs start,end
                 $jxn, $clp1, $clp2, $microH, $insSeq, $jxnSeq, $refRegs); 
    }
  
    # if sufficient non-matching read pairs remain, try to find another junction
    if(@nonJxnKs >= $minSetPairs){
        if (@ks == @nonJxnKs) { # prevent infinite loop in odd case that reconstructed junction matched no read pairs
            # splits and reconstructed clips cannot end here
            # at least one k will match and be put into @jxnKs
            print STDERR "$utility find: processJunction ended with no k match: $jxn\n";
            $jxnOnly or processP2Group();
        } else {  
            @ks = @nonJxnKs;
            findJunction();    
        }
    }
}

# handle sets of read pairs without definitive junction sequence
sub processNoJunction {
    my ($ks) = @_;
    $nGaps++;

    # determine the inferred junction type
    my $jxnType;
    if ($chr1 ne $chr2) {
        $jxnType = "trans";
    } elsif($str1 == $str2){
        $jxnType = $str1 ? "invR" : "invF";
    } else {
        my $k = $$ks[0];
        $jxnType = $pairs[$k][3] - $pairs[$k][0] > $libKeys{$pairs[$k][$#{$pairs[$k]}-1]}[1] ?
                   "del" : "dup";
    }  
    
    # collect the evidence by sample (not library)
    my $maxEvid = 0;
    initEvid();
    foreach my $k(@$ks){
        $evid{all}[0]++;
        $evid{all}[4]++;
        $evid{$pairs[$k][$#{$pairs[$k]}-2]}[0]++;
        $evid{$pairs[$k][$#{$pairs[$k]}-2]}[4]++;
        $maxEvid >= $evid{$pairs[$k][$#{$pairs[$k]}-2]}[0] or
            $maxEvid = $evid{$pairs[$k][$#{$pairs[$k]}-2]}[0];
    }   
    $maxEvid >= $minSetPairs and # at least one sample must have enough evidence
        commitSV($ks, "gap", $jxnType);
}

# set up the null evidence for the input samples for the current candidate SV
sub initEvid {
    foreach my $sample ("all", @samples){
        $evid{$sample} = [(0) x 7];
    }
}

sub commitSV{
    my ($svKs, $jxnSource, $jxnType,                                    # information available from all sources
        $jxn, $clp1, $clp2, $microH, $insSeq, $jxnSeq, $refRegs) = @_;  # information avaiable from sequenced junctions
    
    # suppress duplicate sets (shouldn't happen anymore)
    my $svSig = join(":", sort @$svKs);
    if($svSigs{$svSig}){ $nDupSets++; return }
    $svSigs{$svSig}++;
    $nSets++;
        
    # determine and enforce the sample uniqueness of the SV
    my ($nPosSmp, $nCllSmp) = (0, 0);
    foreach my $sample(@samples){
        $evid{$sample}[0] and $nPosSmp++;
        $evid{$sample}[0] >= $minSetPairs and $nCllSmp++;
    }
    $nPosSmp == 1 and $nUniqueStrict++;
    $nCllSmp == 1 and $nUnique++;
    if($uniqueStrict){
        $nPosSmp > 1 and return;
    } elsif($unique){
        $nCllSmp > 1 and return;
    }

    # assemble the evidence
    my @evidence;
    foreach my $sample("all", @samples){
        push @evidence, "$sample:".join(",", @{$evid{$sample}});
    }    
    
    # describe the set, including quality metrics
    my ($minP1, $maxP1, $minP2, $maxP2, $svSze,
        $fracOv, $fracOv1, $fracOv2, $pSet, $pOutlier) = describeSet($svKs, $jxnType);
    
    # print the variant
    $svID++; # NOTE: the name given is not unique, will be repeated over chromosomes
    print $outH join("\t",  $chr1, $minP1-1, $maxP1, "$group\_$svID\_1", $svSze, "+",
                            $chr2, $minP2-1, $maxP2,
                            $jxnSource, $jxnType,
                            $svSig,
                            $nPosSmp, $nCllSmp,
                            $fracOv, $fracOv1, $fracOv2, $pSet, $pOutlier,
                            @evidence), "\n";
    print $outH join("\t",  $chr2, $minP2-1, $maxP2, "$group\_$svID\_2", $svSze, "+",
                            $chr1, $minP1-1, $maxP1,
                            $jxnSource, $jxnType,
                            $svSig,
                            $nPosSmp, $nCllSmp,
                            $fracOv, $fracOv1, $fracOv2, $pSet, $pOutlier,
                            @evidence), "\n";
    
    # print an index file of the anomalous fragments that reside in the committed sets
    foreach my $k(@$svKs){
        my $maxI = $#{$pairs[$k]};
        my $pairType = $pairs[$k][$maxI];
        $libKey = $pairs[$k][$maxI-1];
        &{$idxSubs{$pairType}}($chr1, $str1, $chr2, $str2, @{$pairs[$k]});
    }  
}

# describe the set for limits, sizes and quality parameters
sub describeSet{  
    my ($svKs, $jxnType) = @_;
    
    # determine the limits of the mapped position of the reads on each side of the junction
    # calculate the structural variant (SV) event size
    my ($minP1, $maxP1, $minP2, $maxP2, $sumSvSze) = (1e9, 0, 1e9, 0, 0);
    foreach my $k(@$svKs){
        $minP1 <= $pairs[$k][0] or $minP1 = $pairs[$k][0];
        $maxP1 >= $pairs[$k][0] or $maxP1 = $pairs[$k][0];
        $minP2 <= $pairs[$k][3] or $minP2 = $pairs[$k][3];
        $maxP2 >= $pairs[$k][3] or $maxP2 = $pairs[$k][3];
        if ($jxnType eq "del") {
            $sumSvSze += $libKeys{$pairs[$k][$#{$pairs[$k]}-1]}[1] - ($pairs[$k][3] - $pairs[$k][0]);
        } elsif ($jxnType eq "dup") {
            $sumSvSze += $libKeys{$pairs[$k][$#{$pairs[$k]}-1]}[1] + ($pairs[$k][0] - $pairs[$k][3]);
        } elsif ($jxnType =~ m/inv/) {
            $sumSvSze += ($pairs[$k][3] - $pairs[$k][0]);
        }
    }
    my $svSze = int($sumSvSze / @$svKs + 0.5);
    
    # calculate extension and overlap region estimates for each fragment in set
    my ($sumOv, %ext) = (0);
    foreach my $k(@$svKs){
        $ext{$k} = ($str1 ? $pairs[$k][0] - $minP1 : $maxP1 - $pairs[$k][0]) + 
                   ($str2 ? $pairs[$k][3] - $minP2 : $maxP2 - $pairs[$k][3]); 
        $sumOv += $libKeys{$pairs[$k][$#{$pairs[$k]}-1]}[1] - $ext{$k};
    }
    
    # calculate overlap estimates for entire set
    my $setOv  = int($sumOv / @$svKs + 0.5);
    #my $setOv_ = $setOv < 1 ? 1 : $setOv;     
    my $fracOv = ($setOv  + ($maxP1 - $minP1) + ($maxP2 - $minP2)) ?
                  int($setOv  / ($setOv  + ($maxP1 - $minP1) + ($maxP2 - $minP2)) * 100) / 100 : 1; 
    #my $fracOv1   = $setOv + ($maxP1 - $minP1) ? int($setOv / ($setOv + ($maxP1 - $minP1) * 2) * 100) / 100 : 1; 
    #my $fracOv2   = $setOv + ($maxP2 - $minP2) ? int($setOv / ($setOv + ($maxP2 - $minP2) * 2) * 100) / 100 : 1;
    my $fracOv1   = ($setOv + ($maxP1 - $minP1) * 2) ? int($setOv / ($setOv + ($maxP1 - $minP1) * 2) * 100) / 100 : 1; 
    my $fracOv2   = ($setOv + ($maxP2 - $minP2) * 2) ? int($setOv / ($setOv + ($maxP2 - $minP2) * 2) * 100) / 100 : 1;
    

    # establish set quality metrics using revised fragment size distribution
    my $pProd = 1;
    my $nOutliers = 0; # fragments with p <0.05 tLen
    foreach my $k(@$svKs){
        my $pTLen = $libKeys{$pairs[$k][$#{$pairs[$k]}-1]}[4];
        my $p = $$pTLen{int(($setOv + $ext{$k})/$pTLenBin + 0.5) * $pTLenBin} || $$pTLen{extreme};
        $pProd *= $p;
        $p < 0.05 and $nOutliers++;
    }
    my $pSet = $pProd ** (1/@$svKs); # the "average" probability for observed fragment sizes
    my $pOutlier = (1 - 0.9**@$svKs) ** $nOutliers;
    return ($minP1, $maxP1, $minP2, $maxP2, $svSze,
            $fracOv, $fracOv1, $fracOv2,
            $pSet, $pOutlier);
}


# CODE BELOW WAS SVTOOLS MIMIC OF VAMP PRIOR TO NEW parsePMatrix approach above

## thread position 1 in ascending order, breaking P1s against each other
#sub threadP1xxx {
#    @js = ();
#    my $prevPair;    
#    foreach my $i(sort {$pairs[$a][0] <=> $pairs[$b][0]} 0..$#pairs){
#        # check for break between current rightmost and next fragment
#        if ($prevPair) {
#            if (($pairs[$i][0] - $$prevPair[0]) > ($str1 ?
#                                                   $pairs[$i][1] - $$prevPair[2] :
#                                                   $$prevPair[1] - $pairs[$i][2])){
#                # prevent premature break by checking current rightmost against next+1 fragment
#                if (!$pairs[$i+1] or
#                    ($pairs[$i+1][0] - $$prevPair[0]) > ($str1 ?
#                                                         $pairs[$i+1][1] - $$prevPair[2] :
#                                                         $$prevPair[1] - $pairs[$i+1][2])){
#                    @js >= $minSetPairs and threadP2();
#                    @js = ();
#                    $prevPair = $pairs[$i];   
#                } else {
#                    $prevPair = $pairs[$i+1]; # pair will never break against itself on next $i
#                }
#            } else {
#                $prevPair = $pairs[$i];
#            }
#        } else {
#            $prevPair = $pairs[$i];
#        }
#        push @js, $i;
#    }
#    @js >= $minSetPairs and threadP2();
#}
#
## thread position 2 in ascending order, breaking P2s against each other
#sub threadP2 {
#    %svSigs = (); # track sets to ensure uniqueness of all committed sets    
#    @ks = ();
#    my $prevPair;    
#    foreach my $j(sort {$pairs[$a][3] <=> $pairs[$b][3]} @js){
#        # check for break between current rightmost and next fragment
#        if ($prevPair) {
#            
#                            #$$vs[$i1]{outPos}, $$vs[$i1]{maxDst}, $$vs[$i1]{minDst},
#                            #$$vs[$i2]{outPos}, $$vs[$i2]{maxDst}, $$vs[$i2]{minDst},    
#            if (($pairs[$j][3] - $$prevPair[3]) > ($str2 ?
#                                                   $pairs[$j][4] - $$prevPair[5] :
#                                                   $$prevPair[4] - $pairs[$j][5])){
#                # prevent premature break by checking current rightmost against next+1 fragment
#                if (!$pairs[$j+1] or
#                    ($pairs[$j+1][3] - $$prevPair[3]) > ($str2 ?
#                                                         $pairs[$j+1][4] - $$prevPair[5] :
#                                                         $$prevPair[4] - $pairs[$j+1][5])){
#                    @ks >= $minSetPairs and findJunction();
#                    @ks = ();
#                    $prevPair = $pairs[$j];
#                } else {
#                    $prevPair = $pairs[$j+1]; # pair will never break against itself on next $j
#                }
#            } else {
#                $prevPair = $pairs[$j];
#            }
#        } else {
#            $prevPair = $pairs[$j];
#        }
#        push @ks, $j;
#    }
#    @ks >= $minSetPairs and findJunction();
#}

#
## break P1/P2 groups into non-overlapping as-mapped groups of fragments, using inferred P1s
#sub processP2Group { # after VAMP
#    
#    
#print STDERR "\n";
#my $continue = 0;
#foreach my $k(@ks){
#    $pairs[$k][3] == 30072602 and $continue = 1;
#}
#$continue or return;
#foreach my $k(@ks){
#    print STDERR "processP2Group: $pairs[$k][0], $pairs[$k][3]\n";
#}
##    30072602
#
#    # perform initial break into as-mapped overlap groups
#    my @ovlpKs;    
#    if ($chr1 eq $chr2){
#        my (@minP1s, @maxP2s); 
#        push @{$ovlpKs[0]}, $ks[0];      
#        @minP1s = ($pairs[$ks[0]][0]);
#        @maxP2s = ($pairs[$ks[0]][3]);
#        foreach my $k(@ks[1..$#ks]){ # fragments supplied by threadP2 in P2 sort order
#            
#            if ($pairs[$k][0] > $maxP2s[$#maxP2s]){
#                
#                push @minP1s, 1E9;     
#                push @maxP2s, 0;
#            }
#            push @{$ovlpKs[$#maxP2s]}, $k;
#            $minP1s[$#minP1s] <= $pairs[$k][0] or $minP1s[$#minP1s] = $pairs[$k][0];
#            $maxP2s[$#maxP2s] = $pairs[$k][3];                                
#        }
#
##processP2Group: 30298241, 30072602
##processP2Group: 30299319, 30074135
##processP2Group: 30298027, 30074659
##processP2Group: 30298070, 30074693
##processP2Group: 30299811, 30076651
#      
#print STDERR "minP1s: ".join("\t", @minP1s), "\n";
#print STDERR "maxP2s: ".join("\t", @maxP2s), "\n";
#
#return;
#
#        
#        # ensure that overlap groups themselves do not overlap (possible given lack of P1 ordering)
#        my $oKI = 1; 
#        while ($ovlpKs[$oKI]){
#            if($minP1s[$oKI] <= $maxP2s[$oKI - 1]){
#                push @{$ovlpKs[$oKI - 1]}, @{$ovlpKs[$oKI]}; # merge overlapping overlap groups          
#                splice(@ovlpKs, $oKI, 1); # removes array element, shifts all remaining elements down one      
#            } else {
#                $oKI++;
#            } 
#        }
#        
#    # as-mapped fragment overlap has no meaning for DiffChrom Sets    
#    } else { 
#        @ovlpKs = \@ks;
#    }
#
#    # finally check overlap groups for self-consistency, breaking into smaller groups if needed
#    foreach my $ks(@ovlpKs){
#        @$ks >= $minSetPairs or next;
#        my $p1Sets = checkSetPositions($ks, 0, $str1);
#        foreach my $p1Set(@$p1Sets){
#            @$p1Set >= $minSetPairs or next;
#            my $p2Sets = checkSetPositions($p1Set, 3, $str2);
#            foreach my $p2Set(@$p2Sets){
#                @$p2Set >= $minSetPairs or next;
#                processNoJunction($p2Set);
#            }
#        }
#    }    
#}
#
## check if the fragments in a group are legitimately contained within a single set
#sub checkSetPositions{ 
#    my ($ks, $posI, $str) = @_;
#    my @sets;     
#    my @ks = sort {$pairs[$a][$posI] <=> $pairs[$b][$posI]} @$ks;
#    if(($pairs[$ks[$#ks]][$posI] - $pairs[$ks[0]][$posI]) >
#           ($str ?
#            $pairs[$ks[$#ks]][$posI+1] - $pairs[$ks[0]][$posI+2] :
#            $pairs[$ks[0]][$posI+1] - $pairs[$ks[$#ks]][$posI+2])){
#        breakSet(\@sets, $ks, $posI, $str);
#    } else {
#        push @sets, $ks;
#    }        
#    return \@sets;
#}
#
## if P1 or P2's are too widely separated, break into multiple, potentially overlapping, sets
#sub breakSet{ 
#    my ($sets, $ks, $posI, $str) = @_;
#    my $maxI = @$ks - 1; 
#    my $iOffset = $minSetPairs - 1;
#    my %goodSpans;
#    
#    # find all non-breaking fragment groups within the collection of fragments
#    foreach my $lowI(0..($maxI - $iOffset)){
#        my $lastI;
#        foreach my $highI(($lowI + $iOffset)..$maxI){ 
#            if(($pairs[$ks[$highI]][$posI] - $pairs[$ks[$lowI]][$posI]) >
#                   ($str ?
#                    $pairs[$ks[$highI]][$posI+1] - $pairs[$ks[$lowI]][$posI+2] :
#                    $pairs[$ks[$lowI]][$posI+1] - $pairs[$ks[$highI]][$posI+2])){
#                last;                
#            } else {
#                $lastI = $highI;
#            }
#        }
#        defined $lastI and $goodSpans{"$lowI:$lastI"}++;
#    }
#    
#    # mask non-breaking fragment groups entirely subsumed by other fragment groups
#    my @goodSpans = keys %goodSpans;
#    foreach my $testSpan(keys %goodSpans){
#        my $keepTestSpan = 1;
#        my ($testLow, $testHigh) = split(":", $testSpan);
#        foreach my $refSpan(@goodSpans){
#            my ($refLow, $refHigh) = split(":", $refSpan);
#            if(($testLow > $refLow and $testHigh <= $refHigh) or ($testLow >= $refLow and $testHigh < $refHigh)) {
#                $keepTestSpan = 0;
#                last;
#            }
#        }
#        $goodSpans{$testSpan} = $keepTestSpan;
#    }
#    
#    # save the final distinct groups
#    foreach my $span(@goodSpans){
#        $goodSpans{$span} or next;
#        my ($lowI, $highI) = split(":", $span); 
#        push @$sets, [@$ks[$lowI..$highI]]; 
#    }
#}

1;
