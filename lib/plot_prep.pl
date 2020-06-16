use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $outH 
            $inputDir $tmpDir $outputDir $plotDir
            %scnRegCols %cnvCols %cndFindCols $nCndFindLeadCol);
my ($group, $inputType, $idCol, $typeCol, $sizeCol, $fracCol, $smpCol, @extraCols, @extraLbls);
my $lbls = '#'.join("\t", qw(chrom start end region unused unused
                             group id cnvType cnvSize cnvFrac sample));

# manage options
sub setOptions_plot_prep {
    setOptionValue(\$group,         'group');  
    setOptionValue(\$inputType,     'input-type');    
    setOptionValue(\$idCol,         'id-col');    
    setOptionValue(\$typeCol,       'type-col');
    setOptionValue(\$sizeCol,       'size-col');
    setOptionValue(\$fracCol,       'frac-col');
    setOptionValue(\$smpCol,        'sample-col'); 
    #-----------------------------------------------
    $idCol   and $idCol--;
    $typeCol and $typeCol--;
    $sizeCol and $sizeCol--;
    $fracCol and $fracCol--;
    $smpCol  and $smpCol--;
    if ($inputType) { # just use some other column if cnvFrac is not defined
        if ($inputType eq 'scan') {
            ($idCol, $typeCol, $sizeCol, $fracCol, $smpCol) =
                map { $scnRegCols{$_} } qw(regId cnvType cnvSize cnvFrac idxSmp);
            @extraLbls = qw(idxBin);
            @extraCols = map { $scnRegCols{$_} } @extraLbls;
        } elsif($inputType eq 'depth') {
            $cnvCols{othSample} = $cnvCols{strand};
            ($idCol, $typeCol, $sizeCol, $fracCol, $smpCol) =
                map { $cnvCols{$_} } qw(cnvId cnvType cnvSize cnvSize sample);
            @extraLbls = qw(othSample);
            @extraCols = map { $cnvCols{$_} } @extraLbls;  
        } elsif($inputType eq 'pairs') {
            $cndFindCols{sample} = $cndFindCols{strand};
            $cndFindCols{othSample} = $cndFindCols{jxnSource};
            ($idCol, $typeCol, $sizeCol, $fracCol, $smpCol) =
                map { $cndFindCols{$_} } qw(svID jxnType svSze svSze sample);
            @extraLbls = qw(othSample);
            @extraCols = map { $cndFindCols{$_} } @extraLbls;
        } else {
            die "unrecognized input type: $inputType\n";
        }
    }
    $idCol    =~ m|\D| and die "$error: option --id-col must be an integer number\n";
    $typeCol  =~ m|\D| and die "$error: option --type-col must be an integer number\n";
    $sizeCol  =~ m|\D| and die "$error: option --size-col must be an integer number\n";
    $fracCol  =~ m|\D| and die "$error: option --frac-col must be an integer number\n";    
    $smpCol   =~ m|\D| and die "$error: option --sample-col must be an integer number\n";
    ($idCol and $typeCol and $smpCol and $sizeCol and $fracCol) or
        die "'input-type' or 'id-col,type-col,size-col,frac-col,sample-col' must be provided\n";
}

# main execution block
sub svtools_plot_prep {

    # initialize
    setOptions_plot_prep();    
    print STDERR "$utility plot_prep: " . getTime(), "\n";

    # collect bin data for each collapse file on the input stream
    openInputStream (undef, \$inH);
    openOutputStream(undef, \$outH);
    print $outH join("\t", $lbls, @extraLbls), "\n";
    my @cols = ($idCol, $typeCol, $sizeCol, $fracCol, $smpCol, @extraCols);
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        my $region = "$f[0]:$f[1]-$f[2]";
        if ($inputType and $inputType eq 'pairs') {
            my %samples;
            foreach my $evid(@f[$nCndFindLeadCol+1..$#f]){ # skip "all"
                $evid =~ m/(.+):(\d+?),/;
                my ($sample, $nFrags) = ($1, $2);
                push @{$samples{$nFrags}}, $sample;
            }
            my @nFrags = sort { $b <=> $a } keys %samples;
            $f[$cndFindCols{sample}] = $samples{$nFrags[0]}[0];
            if (@{$samples{$nFrags[0]}} > 1) {
                $f[$cndFindCols{othSample}] = $samples{$nFrags[0]}[1];
            } elsif($nFrags[1] > 0) {
                $f[$cndFindCols{othSample}] = $samples{$nFrags[1]}[0];
            } else {
                $f[$cndFindCols{othSample}] = 'NA';
            }
        } elsif($inputType and $inputType eq 'depth'){
            $f[$cnvCols{othSample}] = 'NA';
        }
        print $outH join("\t", @f[0..2], $region, 0, '.', $group, @f[@cols], ), "\n";
    }
    closeHandles($inH, $outH);
}

1;
