use strict;
use warnings;

#-----------------------------------------------------------------
# define variables
#-----------------------------------------------------------------
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputFile $region $ouputFile
            $inputDir $tmpDir $outputDir $plotDir
            %findCols $nFindLeadCol %cnvCols $nCnvLeadCol %scnRegCols %biasRegCols);
my ($command, $dataType,
    $humanReadable, $countOnly, $minSize, $maxSize, $minCall, $maxCall, $reqSample, 
    $jxnTypes, $minSetPairs, $maxSetPairs, $minPRatio,
    $jxnOnly, $minPos, $maxPos, $condensed,
    $cnTypes, $minBins, $maxUniq, $minDCN,
    $biasTypes, $minSnps, $minReads, $minBias, $maxOtherBias,
    $scnTypes, $minDltSmp, $maxDltOther, $noGapped, $maxOtherGap, $maxLowQual);
my (%jxnTypes, %cnTypes, %biasTypes, %scnTypes);
my ($gZip, $stepsAfter, $stepsBefore);
my ($count, %count, $filtering, $repeating) = (0);
my (%fCol, $maxFCol);

#-----------------------------------------------------------------
# manage options
#-----------------------------------------------------------------
sub setOptions_view {
    setOptionValue(\$command,       'command');
    setOptionValue(\$dataType,      'data-type');
    #-----------------------------------------------
    setOptionValue(\$ouputFile,     'output-file');      
    setOptionValue(\$humanReadable, 'human');
    setOptionValue(\$countOnly,     'count');
    #-----------------------------------------------
    setOptionValue(\$minSize,       'min-size', 0);
    setOptionValue(\$maxSize,       'max-size', 0);
    setOptionValue(\$minCall,       'min-call', 0);
    setOptionValue(\$maxCall,       'max-call', 0);
    setOptionValue(\$minPos,        'min-pos', 0);
    setOptionValue(\$maxPos,        'max-pos', 0);    
    setOptionValue(\$reqSample,     'sample',   '');
    setOptionValue(\$minBins,       'min-bins', 0);    
    setOptionValue(\$minSetPairs,   'min-set-pairs', 0);
    setOptionValue(\$maxSetPairs,   'max-set-pairs', 0);
    setOptionValue(\$minPRatio,     'min-PRatio', 0);    
    #-----------------------------------------------    
    setOptionValue(\$jxnTypes,      'jxn-types', '');  
    setOptionValue(\$jxnOnly,       'jxn-only');
    #setOptionValue(\$unique,        'unique');
    #setOptionValue(\$uniqueStrict,  'unique-strict');
    setOptionValue(\$condensed,       'condensed');
    #----------------------------------------------- 
    setOptionValue(\$cnTypes,       'cn-types', '');
    setOptionValue(\$maxUniq,       'max-uniq', 0);
    setOptionValue(\$minDCN,        'min-delta-CN', 0);
    #----------------------------------------------- 
    setOptionValue(\$biasTypes,     'bias-types', '');
    setOptionValue(\$minSnps,       'min-snps', '');
    setOptionValue(\$minReads,      'min-reads', '');
    setOptionValue(\$minBias,       'min-bias', '');
    setOptionValue(\$maxOtherBias,  'max-other-bias', '');
    #----------------------------------------------- 
    setOptionValue(\$scnTypes,      'scan-types', '');  
    setOptionValue(\$minDltSmp,     'min-dlt-sample');
    setOptionValue(\$maxDltOther,   'max-dlt-other');
    setOptionValue(\$noGapped,      'suppress-gapped');
    setOptionValue(\$maxOtherGap,   'max-other-gap');
    setOptionValue(\$maxLowQual,    'max-low-qual');
    #-----------------------------------------
    $minSize =~ m|\D| and die "$error: min-size must be an integer number of bp\n";     
    $maxSize =~ m|\D| and die "$error: max-size must be an integer number of bp\n";    
    $minCall =~ m|\D| and die "$error: min-call must be an integer number\n";     
    $maxCall =~ m|\D| and die "$error: max-call must be an integer number\n";
    #-----------------------------------------
    checkTypes('junction', $jxnTypes, \%jxnTypes, qw(del dup invF invR trans));
    $minSetPairs =~ m|\D| and die "$error: min-set-pairs must be an integer number\n";    
    $maxSetPairs =~ m|\D| and die "$error: max-set-pairs must be an integer number\n";  
    $minPos =~ m|\D| and die "$error: min-pos must be an integer number\n";     
    $maxPos =~ m|\D| and die "$error: max-pos must be an integer number\n";
    #-----------------------------------------
    checkTypes('cn', $cnTypes, \%cnTypes, qw(loss gain group));
    $minBins =~ m|\D| and die "$error: min-bins must be an integer number of bins\n";
    #-----------------------------------------
    checkTypes('bias', $biasTypes, \%biasTypes, qw(refDup altDup refDel altDel));
    $minSnps =~ m|\D| and die "$error: min-snps must be an integer number of bins\n";
    $minReads =~ m|\D| and die "$error: min-reads must be an integer number of bins\n";   
    #-----------------------------------------
    checkTypes('scn', $scnTypes, \%scnTypes, qw(del dup ins gain loss));
    #-----------------------------------------
    ($inputFile or ($command and $dataType)) or
        die "$error: command and data-type must be specified when data are provided on stdin\n";
}
sub checkTypes {
    my ($typeType, $types, $hash, @types) = @_;
    if ($types) {
        %$hash = map {$_ => 1} split(",", $types);
        my %known = map { $_ => 1 } @types;
        foreach my $type(keys %$hash){
            $known{$type} or die "$error: unknown $typeType type: $type\n";
        }
    }    
}

#-----------------------------------------------------------------
# main execution block
#-----------------------------------------------------------------
sub svtools_view {
    ($inputFile, $region) = @_;
    setOptions_view();

    # open IO streams
    if($inputFile){
        $inputFile =~ m/svtools\.(\w+)\.(\w+)\..+\.bgz/ or die "$error: could not parse file name: $inputFile\n";
        ($command, $dataType) = ($1, $2);
        if ($region) {
            openInputStream(undef, \$inH, undef, "tabix $inputFile $region");
        } else {
            openInputStream(undef, \$inH, undef, "cat $inputFile | gunzip -c");
        }
    }  else {
        openInputStream(undef, \$inH);
    }
    openOutputStream($ouputFile, \$outH, undef, undef);
    
    # set the type, and determine filtering
    my $viewSub;
    $repeating = !($humanReadable or $countOnly);
    if (($command eq 'cat' or $command eq 'find') and $dataType eq 'sets') {
        $viewSub = \&viewSets;
        $filtering = ($minSize or $maxSize or $minCall or $maxCall or $reqSample or
                      $jxnTypes or $minSetPairs or $maxSetPairs or
                      $jxnOnly or $minPos or $maxPos);
        %fCol = %findCols;
        $maxFCol = $nFindLeadCol;
    } elsif ($command eq 'call' and $dataType eq 'cnvs'){
        $viewSub = \&viewCnvs;
        $filtering = ($minSize or $maxSize or $minCall or $maxCall or $reqSample or
                      $cnTypes or $minBins or $maxUniq or $minDCN or $minPRatio or
                      $minSetPairs or $maxSetPairs);
        %fCol = %cnvCols;
        $maxFCol = $nCnvLeadCol;
    } elsif ($command eq 'collapse_bias' and $dataType eq 'regions'){
        $viewSub = \&viewBias;
        $filtering = ($minSize or $maxSize or $minCall or $maxCall or $minPos or $maxPos or $reqSample or
                      $biasTypes or $minSnps or $minReads or $minBias or $maxOtherBias or $minPRatio);
        %fCol = %biasRegCols;          
    } elsif ($command eq 'collapse' and $dataType eq 'regions'){
        $viewSub = \&viewScan;
        $filtering = ($minSize or $maxSize or $minCall or $maxCall or $minPos or $maxPos or $reqSample or
                      $minBins or
                      $scnTypes or $minDltSmp or $maxDltOther or
                      $noGapped or $maxOtherGap or $maxLowQual);
        %fCol = %scnRegCols;
    } else {
        die "$error: unrecognized data-type: $command $dataType\n";
    }
    
    # display the results    
    while (my $line = <$inH>) {
        if ($repeating and !$filtering) {
            print $line;
        } else {
            chomp $line;
            &$viewSub($line);
        }
    }
    $countOnly and print "$count rows".($dataType eq 'sets' ? ", ".scalar(keys %count)." junctions\n" : "\n");
    closeHandles($inH, $outH);
}

#-----------------------------------------------------------------
# parse the output from find...cat
#-----------------------------------------------------------------
sub viewSets {
    my ($line) = @_;
    my @f = split("\t", $line);
    if ($filtering) {
        if($minSize){
            $f[$fCol{chrom1}] ne $f[$fCol{chrom2}] and return;
            abs($f[$fCol{svSze}]) < $minSize and return;
        }
        if($maxSize){
            $f[$fCol{chrom1}] ne $f[$fCol{chrom2}] and return;
            abs($f[$fCol{svSze}]) > $maxSize and return;
        }        
        $minCall and ($f[$fCol{nCllSmp}] >= $minCall or return);
        $maxCall and ($f[$fCol{nCllSmp}] <= $maxCall or return);
        $jxnTypes and ($jxnTypes{$f[$fCol{jxnType}]} or return);        
        if($minSetPairs or $maxSetPairs){
            my @frags = split(":", $f[$fCol{frags}]);
            my $nFrags = scalar(@frags);
            $minSetPairs and ($nFrags >= $minSetPairs or return);
            $maxSetPairs and ($nFrags <= $maxSetPairs or return);
        }
        $minPos and ($f[$fCol{nPosSmp}] >= $minPos or return);
        $maxPos and ($f[$fCol{nPosSmp}] <= $maxPos or return);
        #$uniqueStrict and ($f[$fCol{nPosSmp}] == 1 or return);
        #$unique and ($f[$fCol{nCllSmp}] == 1 or return);
        #$reqSample or $jxnOnly
    }            
#chrom1 minPThis maxPThis svID svSze strand
#chrom2 minPThat maxPThat
#jxnSource jxnType
#frags
#nPosSmp nCllSmp
#fracOvlp fracOvlp1 fracOvlp2
#pSet pOutlier
    
    if ($countOnly) {
        $count++;
        $f[$fCol{svID}] =~ m/(.+)_\d+/ and $count{$1}++;
    } elsif($humanReadable) {
        print "-" x 80, "\n";
        print join("\t", $f[$fCol{svID}], commify($f[$fCol{svSze}])."bp", $f[$fCol{jxnType}]), "\n";
        print "$f[$fCol{chrom1}]:".commify($f[$fCol{minPThis}])."-".commify($f[$fCol{maxPThis}]).
              " with ".
              "$f[$fCol{chrom2}]:".commify($f[$fCol{minPThat}])."-".commify($f[$fCol{maxPThat}])."\n";
        if ($f[$fCol{chrom1}] eq $f[$fCol{chrom2}]) {
            my $start = min($f[$fCol{minPThis}], $f[$fCol{minPThat}]);
            my $end   = max($f[$fCol{maxPThis}], $f[$fCol{maxPThat}]);                    
            print "$f[$fCol{chrom1}]:$start-$end\n";
        }
        print "FRAGS\t$f[$fCol{frags}]\n";
        printPaddedMatrix(
            [qw(nPosSmp nCllSmp fracOvlp fracOvlp1 fracOvlp2 pSet pOutlier)],
            [@f[$fCol{nPosSmp}..$fCol{pOutlier}]]              
        );
        printPaddedMatrix(
            #[qw(sample evidence)],
            [split(/:|,/, $f[$maxFCol])],
            ["--------", ("-") x 7],
            map { [split(/:|,/, $f[$_])] } ($maxFCol+1)..$#f
        );
    } elsif($condensed) {
        $f[$fCol{chrom1}] eq $f[$fCol{chrom2}] or return;        
        $f[$fCol{svID}] =~ m/(.+)_(\d)/ or die "bad svID in input\n";
        $2 == 1 or return;
        my $svId = $1;
        my $start = min($f[$fCol{minPThis}], $f[$fCol{minPThat}]);
        my $end   = max($f[$fCol{maxPThis}], $f[$fCol{maxPThat}]);
        print join("\t", 
            $f[$fCol{chrom1}], $start, $end, $svId,
            @f[$fCol{svSze}..$fCol{strand}],
            @f[$fCol{jxnSource}..$#f]
        ), "\n";
    } else {
        print $line, "\n";
    }  
}

#-----------------------------------------------------------------
# parse the output from bin...call
#-----------------------------------------------------------------
sub viewCnvs {
    my ($line) = @_;
    my @f = split("\t", $line);  
    if ($filtering) {
        if($minSize){
            $f[$fCol{end}] - $f[$fCol{start}] < $minSize and return;
        }
        if($maxSize){
            $f[$fCol{end}] - $f[$fCol{start}] > $maxSize and return;
        }
        $minCall and ($f[$fCol{maxNSamples}] >= $minCall or return);
        $maxCall and ($f[$fCol{maxNSamples}] <= $maxCall or return);
        $reqSample and ($f[$fCol{sample}] eq $reqSample or return);
        if($cnTypes){
            my $cnT = $f[$fCol{cnvType}];
            $cnT eq 'cnv' and $cnT = 'group';
            $cnTypes{$cnT} or return;
        }        
        $minBins and ($f[$fCol{nBins}] >= $minBins or return);
        $maxUniq and ($f[$fCol{uniqueness}] <= $maxUniq or return);
        if ($minDCN or $minPRatio) {
            $f[$fCol{cnvType}] eq 'cnv' and return;        
            my $sampleIndex = $f[$fCol{sampleIndex}];
            if ($minDCN) {
                my $dcn = $f[$nCnvLeadCol + $sampleIndex];
                abs($dcn) >= $minDCN or return; # abs(), since neg dCN is a loss
            }
            if ($minPRatio) {
                my $nSamples = (@f - $nCnvLeadCol) / 2;
                my $pRatio = $f[$nCnvLeadCol + $sampleIndex + $nSamples];
                $pRatio >= $minPRatio or return; # not abs(), neg pRatio favors neutral
            }
        }   
    }            
# chrom start end cnvId nBins strand
# cnvType cnvSize sample sampleIndex
# uniqueness maxNSamples binsList
#   samples x mean(delta_copy_number_raw)
#   samples x PRatio(cnvType vs. neutral)
    if ($countOnly) {
        $count++;
    } elsif($humanReadable) {
        print "-" x 80, "\n";
        print join("\t", $f[$fCol{sample}], "id:$f[$fCol{cnvId}]", commify($f[$fCol{cnvSize}])." bp", $f[$fCol{cnvType}]), "\n";
        print "$f[$fCol{chrom}]:".commify($f[$fCol{start}])."-".commify($f[$fCol{end}]), "\n";
        printPaddedMatrix(
            [qw(nBins uniqueness maxNSamples)],
            [@f[$fCol{nBins}, $fCol{uniqueness}..$fCol{maxNSamples}]]              
        );
        my $sampleIndex = $f[$fCol{sampleIndex}];
        my $nSamples = (@f - $nCnvLeadCol) / 2;
        my @dcn    = @f[$nCnvLeadCol..$nCnvLeadCol+$nSamples-1];
        my @pRatio = @f[$nCnvLeadCol+$nSamples..$#f];
        if ($sampleIndex >= 0) {
            $dcn[$sampleIndex] = $dcn[$sampleIndex]."*";
            $pRatio[$sampleIndex] = $pRatio[$sampleIndex]."*";
        }
        print "samples...\n";
        printPaddedMatrix( \@dcn, \@pRatio );     
        
    } else {
        print $line, "\n";
    }  
}
#-----------------------------------------------------------------
# parse the output from genotype...collapse_bias
#-----------------------------------------------------------------
sub viewBias {
    my ($line) = @_;
    my @f = split("\t", $line);                
    if ($filtering) {
        if($minSize){ $f[$fCol{cnvSize}] < $minSize and return }        
        if($maxSize){ $f[$fCol{cnvSize}] > $maxSize and return }
        $minCall and ($f[$fCol{nCllSmp}] >= $minCall or return);
        $maxCall and ($f[$fCol{nCllSmp}] <= $maxCall or return);
        $minPos and  ($f[$fCol{nPosSmp}] >= $minPos or return);
        $maxPos and  ($f[$fCol{nPosSmp}] <= $maxPos or return);        
        $reqSample and ($f[$fCol{idxSmp}] eq $reqSample or return);         
        if($biasTypes){ $biasTypes{$f[$fCol{biasType}]} or return }
        $minSnps and ($f[$fCol{nSnps}] >= $minSnps or return);        
        $minReads and (($f[$fCol{nRef}] + $f[$fCol{nAlt}]) >= $minReads or return);
        $minPRatio and ($f[$fCol{maxPRatioIdx}] >= $minPRatio or return);
        $minBias and (abs($f[$fCol{bias}]) >= $minBias or return);
        $maxOtherBias and (abs($f[$fCol{maxOtherBias}]) <= $maxOtherBias or return);   
    }         
#   chrom start end regId score strand
#	nPosSmp nCllSmp
#	idxSmp maxPRatioIdx maxPRatioOther
#   start_ end_ cnvId nSnps nRef nAlt biasType
#   cnvSize bias maxOtherBias
    if ($countOnly) {
        $count++;
    } elsif($humanReadable) {
        print "-" x 80, "\n";        
        print join("\t", $f[$fCol{idxSmp}], "id:$f[$fCol{regId}]", commify($f[$fCol{cnvSize}])." bp", $f[$fCol{biasType}]), "\n";        
        print "$f[$fCol{chrom}]:".commify($f[$fCol{start}])."-".commify($f[$fCol{end}]), "\n";        
        printPaddedMatrix(
            [qw(nSnps
                nPosSmp nCllSmp
                maxPRatioIdx maxPRatioOther
                nRef nAlt bias maxOtherBias)],
            [@f[$fCol{nSnps},
                $fCol{nPosSmp}, $fCol{nCllSmp},
                $fCol{maxPRatioIdx}, $fCol{maxPRatioOther},
                $fCol{nRef}, $fCol{nAlt}, $fCol{bias}, $fCol{maxOtherBias}]]              
        );
    } else {
        print $line, "\n";
    }  
}

#-----------------------------------------------------------------
# parse the output from scan...collapse
#-----------------------------------------------------------------
sub viewScan {
    my ($line) = @_;
    my @f = split("\t", $line);
    if ($filtering) {
        if($minSize){ abs($f[$fCol{cnvSize}]) < $minSize and return }        
        if($maxSize){ abs($f[$fCol{cnvSize}]) > $maxSize and return }
        $minCall and ($f[$fCol{nCllSmp}] >= $minCall or return);
        $maxCall and ($f[$fCol{nCllSmp}] <= $maxCall or return);
        $minPos and  ($f[$fCol{nPosSmp}] >= $minPos or return);
        $maxPos and  ($f[$fCol{nPosSmp}] <= $maxPos or return);
        $reqSample and ($f[$fCol{idxSmp}] eq $reqSample or return);        
        $minBins and ($f[$fCol{nRegBins}] >= $minBins or return);        
        if($scnTypes){ $scnTypes{$f[$fCol{cnvType}]} or return } 
        $minDltSmp    and ($f[$fCol{maxPosDltIdx}] >= $minDltSmp or return);
        $maxDltOther  and ($f[$fCol{maxPosDltOther}] <= $maxDltOther or return);        
        $minSetPairs  and ($f[$fCol{nFrags}] >= $minSetPairs or return);
        $maxSetPairs  and ($f[$fCol{nFrags}] <= $maxSetPairs or return);  
        $noGapped     and $f[$fCol{isGap}] and return;
        $maxOtherGap  and ($f[$fCol{nOtherGap}] <= $maxOtherGap or return);
        $maxLowQual   and ($f[$fCol{fracLowQual}] <= $maxLowQual or return);
    }         
#   chrom start end regId nRegBins strand 
#	nPosSmp nCllSmp
#	idxSmp maxCnvDltIdx maxPosDltIdx maxCnvDltOther maxPosDltOther
#	isGap nOtherGap
#   idxBin nBins nFrags
#	fracLowQual tLens isGap n_log_p deltaCDF_abs deltaCDF_sgn
#	deltaCDF_alt cnvType cnvSize cnvFrac
    if ($countOnly) {
        $count++;
    } elsif($humanReadable) {
        print "-" x 80, "\n";        
        print join("\t", $f[$fCol{idxSmp}], "id:$f[$fCol{regId}]", commify($f[$fCol{cnvSize}])." bp", $f[$fCol{cnvType}]), "\n";        
        print "$f[$fCol{chrom}]:".commify($f[$fCol{start}])."-".commify($f[$fCol{end}]), "\n";        
        printPaddedMatrix(
            [qw(nFrags cnvFrac nRegBins
                nPosSmp nCllSmp
                maxDltIdx maxDltOther
                fracLowQual isGap nOtherGap)],
            [@f[$fCol{nFrags},  $fCol{cnvFrac},  $fCol{nRegBins},
                $fCol{nPosSmp}, $fCol{nCllSmp},
                $fCol{maxPosDltIdx}, $fCol{maxPosDltOther},
                $fCol{fracLowQual}, $fCol{isGap}, $fCol{nOtherGap}]]              
        );
    } else {
        print $line, "\n";
    }  
}
#-----------------------------------------------------------------
# print subs
#-----------------------------------------------------------------
sub printPaddedMatrix {
    my (@rows) = @_;
    my @maxLen = (0) x 100;
    foreach my $row(@rows){
        foreach my $i(0..$#$row){
            my $len = length($$row[$i]);
            $maxLen[$i] >= $len or $maxLen[$i] = $len;
        }
    }
    foreach my $row(@rows){
        my @col;
        foreach my $i(0..$#$row){
            my $len = length($$row[$i]);
            push @col, (" " x ($maxLen[$i] - $len)).$$row[$i];
        }
        print join("  ", @col), "\n";
    }
}

1;
