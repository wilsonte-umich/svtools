use strict;
use warnings;
use List::Util qw(shuffle);

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $tmpH
            $inputDir $tmpDir $outputDir $plotDir
            %scnPosCols);
my %c = %scnPosCols;
my (@globs, $maxMem,$sample, $pdfDir,
    $minCnvDelta, $minCnvBins,
    $cnvSizeStep, $minCnvSize, $maxCnvSize,
    $region);
my ($minTLen, $maxTLen);
my ($scnH, $cnvH);
my ($prevChrom, $prevBin, @cnvBins);

# manage options
sub setOptions_crunch {
    setOptionValue(\$sample,     'sample');
    setOptionValue(\$inputDir,  'input-dir');
	setOptionValue(\$pdfDir,     'pdf-dir');  
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);
    setOptionValue(\$outputDir,  'output-dir');
	setOptionValue(\$minCnvDelta,'min-cnv-delta', 5);
    setOptionValue(\$minCnvBins, 'min-cnv-bins',  5);
    setOptionValue(\$cnvSizeStep,'cnv-size-step', 500);
    setOptionValue(\$minCnvSize, 'min-cnv-size',  2000);
    setOptionValue(\$maxCnvSize, 'max-cnv-size',  10000);
    setOptionValue(\$region,     'region');
    #-------------------------------------
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
	-d $pdfDir or die "$error: $pdfDir does not exist or is not a directory\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    $minCnvBins =~ m|\D| and die "$error: min-cnv-bins must be an integer number\n";
    $minCnvSize =~ m|\D| and die "$error: min-cnv-size must be an integer number of bp\n";
    $maxCnvSize =~ m|\D| and die "$error: max-cnv-size must be an integer number of bp\n";
    $minCnvSize % $cnvSizeStep and die "$error: min-cnv-size must be a multiple of cnv-size-step\n";
    $maxCnvSize % $cnvSizeStep and die "$error: max-cnv-size must be a multiple of cnv-size-step\n";    
    $region and ($region =~ m/(chr.+):\d+-\d+/ or die "region must be in form chr:start-end");
    #-------------------------------------
    $minTLen = -$maxCnvSize;
    $maxTLen = 2 * $maxCnvSize;
}

# main execution block
sub svtools_crunch {

    # initialize
    setOptions_crunch();    
    print STDERR "$utility crunch: " . getTime(), "\n";

    # find runs
	print STDERR "finding runs of deviant bins\n";    
    my $binFile = getInFile('scan', 'positions', $sample, 'bgz');
    if ($region) {
        openInputStream("tabix $binFile $region", \$scnH);	 
    } else {
        openInputStream("gunzip -c $binFile", \$scnH);	        
    }
    my $tmpFile = getTmpFile("cnvs", $sample, "gz");
    openOutputStream($tmpFile, \$tmpH, $TRUE);
    find_runs();
    closeHandles($scnH, $tmpH);

    # assign CNV values
	print STDERR "applying area minimization to find CNV values\n";
	my $cnvFile = getOutFile("cnvs", $sample, "bgz"); 
	unless($region){
        openOutputStream($cnvFile, \$cnvH, undef, undef,
                "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    }
	describe_cnvs($tmpFile);
	closeHandles($cnvH); 
 
    # finish up
    unlink $tmpFile;
	unless($region){
        system("tabix -p bed $cnvFile");
    }
}

# find runs of deviant bins
sub find_runs {
    my $inRun;
    while (my $line = <$scnH>) {
        chomp $line;
		my @f = split("\t", $line);
        my $isDev = (abs($f[$c{deltaCDF_sgn}]) >= $minCnvDelta);
        $inRun and !$isDev and commitCNV();
        $inRun = $isDev;
        $isDev and push @cnvBins, \@f;
    }
    commitCNV();
}
sub commitCNV {
	@cnvBins or return;
    my $nBins = @cnvBins;
    if($nBins >= $minCnvBins){
        my $midI = int($nBins / 2);
        print $tmpH  join("\t",
            $cnvBins[0][$c{chrom}], $cnvBins[0][$c{start}], $cnvBins[$#cnvBins][$c{end}],
            $sample, $nBins, '.',
            @{$cnvBins[$midI]}[$c{nFrags},$c{fracLowQual}..$c{deltaCDF_sgn}]
        )."\n";
    }
	@cnvBins = ();	
}

# analyze CNVs in R and commit the results
sub describe_cnvs {
	my ($tmpFile) = @_;
    $ENV{sample}      = $sample;
	$ENV{pdfDir} 	  = $pdfDir;
	$ENV{cnvFile}     = $tmpFile;
	$ENV{cnvSizeStep} = $cnvSizeStep;
	$ENV{minCnvSize}  = $minCnvSize;
	$ENV{maxCnvSize}  = $maxCnvSize;
    open(my $rH, "-|", "Rscript $libDir/crunch.R") or die "could not open crunch.R stream: $!\n";
    while (my $line = <$rH>) {
		if ($region) {
			print $line;
		} else {
			print $cnvH $line;
		}
    }
    closeHandles($rH);
}

1;
