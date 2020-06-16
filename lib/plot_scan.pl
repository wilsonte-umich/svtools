use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $tmpH 
            $inputDir $tmpDir $outputDir $plotDir
            %pPrpCols %scnPosCols);
my %p = %scnPosCols;
my %c = %pPrpCols;
$c{idxBin} = (sort {$b <=> $a} values %c)[0] + 1;
my (@globs, $group, $pdfDir, $centerWidth, $widthPadding, $bamGlob, $maxCnvSize, $minQual);
my ($minTLen, $maxTLen) = @_;
my (@files, %samples, @samples);
my ($posH, $datH, $spnH, $clpH, $lstH, $bamH);

# manage options
sub setOptions_plot_scan {
    setOptionValue(\$group,         'group');
    setOptionValue(\$plotDir,       'plot-dir');
    setOptionValue(\$pdfDir,        'pdf-dir');  
    setOptionValue(\$tmpDir,        'tmp-dir',      '/tmp');
    setOptionValue(\$centerWidth,   'center-width', 5);
    setOptionValue(\$bamGlob,       'bam-glob');
    setOptionValue(\$maxCnvSize,    'max-cnv-size', 10000);
    setOptionValue(\$minQual,       'min-qual',     5);
    #------------------------------------------------------
    -d $plotDir or die "$error: $plotDir does not exist or is not a directory\n";
    -d $pdfDir or die "$error: $pdfDir does not exist or is not a directory\n";
    -d $tmpDir  or die "$error: $tmpDir does not exist or is not a directory\n";    
    $centerWidth =~ m|\D| and die "$error: center-width must be an integer number\n";
    $minQual =~ m|\D| and die "$error: min-qual must be an integer MAPQ\n";    
    $outputDir = $plotDir;
    $widthPadding = ($centerWidth - 1) / 2;
    #------------------------------------------------------
    $minTLen = -$maxCnvSize;
    $maxTLen = 2 * $maxCnvSize;
}

# main execution block
sub svtools_plot_scan {
    (@globs) = @_;
    $globs[0] or die "$error: must supply list of 'svtools.scan.positions' files on command line\n";

    # initialize
    setOptions_plot_scan();    
    print STDERR "$utility plot_scan: " . getTime(), "\n";
    foreach my $glob(@globs){ push @files, glob($glob) }
    my ($rand, $n) = (int(rand(1e6)), 1);
    foreach my $file(@files){
        $file =~ m/svtools\.scan\.positions\.(.+)\.bgz$/ or die "not an scan positions file: $file\n";
        $samples{$1}++;
    }
    @samples = keys %samples;

    # collect bin data for each collapse file on the input stream
    print STDERR "collecting position data\n";
    my $plotFile = getTmpFile("plot_scan", 'data', "gz");
    openOutputStream($plotFile, \$datH, $TRUE);
    my $spanFile = getTmpFile("plot_scan", 'spans', "gz");
    openOutputStream($spanFile, \$spnH, $TRUE);
    my $lstFile = getOutFile("list", $group, "txt");
    openOutputStream($lstFile, \$lstH, $FALSE);    
    openInputStream(undef, \$clpH);
    my $header = <$clpH>;
    print $lstH $header;    
    while (my $clp = <$clpH>) {
        print $lstH $clp;
        chomp $clp;
        my @clp = split("\t", $clp);
        my $width = $clp[$c{end}] - $clp[$c{start}] - 1;
        my $start = $clp[$c{start}] - $width * $widthPadding;
        my $end = $clp[$c{end}] + $width * $widthPadding;
        my $plotRegion = "$clp[$c{chrom}]:$start-$end";
        $clp[$c{region}] =~ s/:/-/;
        print STDERR "  $clp[$c{grpId}]\t";
        #print join("\t", @clp), "\n";
        foreach my $file(@files){
            open my $tbxH, "-|", "tabix $file $plotRegion" or
                die "could not open tabix stream: tabix $file $plotRegion\n";
            print STDERR ".";
            while (my $tbx = <$tbxH>) {
                chomp $tbx;
                my @tbx = split("\t", $tbx);
                print $datH join("\t",
                    (map {$clp[$c{$_}]} qw(chrom region group grpId cnvType cnvSize cnvFrac sample idxBin)),                                 
                    (map {$tbx[$p{$_}]} qw(end sample n_log_p deltaCDF_sgn nFrags))), "\n";
            }
            close $tbxH;
            $file =~ m/svtools\.scan\.positions\.(.+)\.bgz$/;
            my $sample = $1;
            addBamSpans($clp[$c{grpId}], $clp[$c{sample}], $sample, $clp[$c{chrom}], $start, $end); 
        }
        print STDERR "\n";
    }
    closeHandles($datH, $spnH, $clpH, $lstH);

    # make the plots
    print STDERR "making plots\n";
    $ENV{plotFile}     = $plotFile;
    $ENV{spanFile}     = $spanFile;
    $ENV{plotDir}      = $plotDir;
    $ENV{pdfDir}       = $pdfDir;
    $ENV{centerWidth}  = $centerWidth;
    $ENV{minTLen}      = $minTLen;
    $ENV{maxTLen}      = $maxTLen;
    $ENV{minQual}      = $minQual;
    system("Rscript $libDir/plot_scan.R");
    
    # finish up
    unlink $plotFile;
    unlink $spanFile;
}

sub addBamSpans {
    my ($grpId, $idxSample, $sample, $chrom, $start, $end) = @_;
    my $bam = $bamGlob;
    $bam =~ s/SAMPLE/$sample/;
    my @files = glob($bam);
    my $bamStart = $start - $maxTLen;
    my $bamEnd = $end + $maxTLen;
    my $bamRegion = "$chrom:$bamStart-$bamEnd";  
    open my $bamH, "-|", "samtools merge -u -R $bamRegion - @files | samtools view -" or
        die "could not open bam stream: $bam $bamRegion\n";
    my %enc;
    while (my $bam = <$bamH>) {
        chomp $bam;
        my @aln = split("\t", $bam, 12);
        $aln[11] =~ m/RG:Z:(\S+)/ or die "missing read group in line:\n$bam";
        my $lib = $1;
        if(($aln[1] & 0x930) == 0x020 and # forward read of FR pair, not secondary or supplemental
            $aln[6] eq '=' and            # on same chromosome		
            $aln[8] >= $minTLen and $aln[8] <= $maxTLen){ # within interrogation size range
                my $fPos = $aln[3];
                my $rPos = $fPos + $aln[8];
                my ($leftPos, $rghtPos) = ($fPos <= $rPos) ? ($fPos, $rPos) : ($rPos, $fPos);
                my $pairKey = "$lib,$fPos,$aln[8]";
                $enc{$pairKey} and next;
                $enc{$pairKey}++;
                print $spnH join("\t", $grpId, $idxSample, $sample, $leftPos, $rghtPos, $aln[8], $aln[4]), "\n";
        }    
    }
    close $bamH;
}

1;
