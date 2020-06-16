use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH
            $inputDir $tmpDir $outputDir);
require "$libDir/exclude_subs.pl";
my (@globs, $cellRangerDir, $sample, $genome, $includeY); 

# manage options
sub setOptions_cluster10X {
    setOptionValue(\$inputDir,  'input-dir');
    setOptionValue(\$cellRangerDir, 'cell-ranger-dir');
    setOptionValue(\$sample,    'sample');
    setOptionValue(\$genome,    'genome');
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$includeY,  'include-Y');
    #-------------------------------------
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    -d $cellRangerDir or die "$error: $cellRangerDir does not exist or is not a directory\n";    
    $outputDir or $outputDir = $inputDir;
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
}

# main execution block
sub svtools_cluster10X {
    (@globs) = @_;
    $globs[0] and die "$error: cluster10X does not take any input (files are read from disk)\n";

    # initialize
    setOptions_cluster10X();    
    print STDERR "$utility cluster10X: " . getTime(), "\n";

    # normalize cluster cell data
    print STDERR "$utility cluster10X: calculating Z scores and performing hierarchical clustering\n";
    $ENV{LIB_DIR}   = $libDir;
    $ENV{CR_DIR}    = $cellRangerDir;
    $ENV{IN_DIR}    = $inputDir;
    $ENV{OUT_DIR}   = $outputDir;
    $ENV{SAMPLE}    = $sample;
    $ENV{GENOME}    = $genome;
    $ENV{INCLUDE_Y} = $includeY ? 'TRUE' : 'FALSE'; 
    system("Rscript $libDir/cluster10X.R") and die "$error: cluster10X.R returned an error\n";  
}

1;
