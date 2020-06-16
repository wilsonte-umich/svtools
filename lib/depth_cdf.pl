use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            %binCols %binColOffsets $nBinLeadCol
            %binFileCols $nbinFileLeadCol);
my (@globs, $group);

# manage options
sub setOptions_depth_cdf {
    setOptionValue(\$group,     'group');
    setOptionValue(\$inputDir,  'input-dir');
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$tmpDir,    'tmp-dir',         '/tmp');
    #-------------------------------------
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
    -d $tmpDir    or die "$error: $tmpDir does not exist or is not a directory\n";
    -d $inputDir  or die "$error: $inputDir does not exist or is not a directory\n";
}

# main execution block
sub svtools_depth_cdf {
    # take no files as input, collected based on group and inputDir

    # initialize
    setOptions_depth_cdf();    
    print STDERR "$utility call: " . getTime(), "\n";
    my $smpFile = getInFile("call", "samples", $group, "txt");
    my $binFile = getInFile("call", "bins",    $group, "bgz");
    my $tmpFile = getTmpFile("depth_cdf", $group, "txt");    
    my $rFile   = getOutFile("depth_cdf", $group, "R");
    my $jpgFile = getOutFile("depth_cdf", $group, "jpg");

    # parse the data for R
    print STDERR "$utility call: plotting bin count CDF\n";
    openOutputStream($tmpFile, \my $tmpH);
    my @samples = getSamples($smpFile);
    print $tmpH join("\t", @samples), "\n";
    openInputStream($binFile, \$inH, $FALSE, $FALSE, $TRUE);
    my @col =  $nbinFileLeadCol +  $binColOffsets{smpCovNrm}      * @samples ..
              ($nbinFileLeadCol + ($binColOffsets{smpCovNrm} + 1) * @samples) - 1;
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        print $tmpH join("\t", @f[@col]), "\n";
    }
    closeHandles($tmpH, $inH);
    
    # make the plot
    $ENV{tmpFile} = $tmpFile;
    $ENV{rFile}   = $rFile;
    $ENV{jpgFile} = $jpgFile;
    system("Rscript $libDir/depth_cdf.R");
    
    # clean up
    unlink $tmpFile;
}

1;
