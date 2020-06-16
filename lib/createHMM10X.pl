use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH
            $inputDir $tmpDir $outputDir $plotDir);
my (@globs, $weightPerCell);

# manage options
sub setOptions_createHMM10X {
    setOptionValue(\$weightPerCell, 'weight-per-cell', 10); 
    #-------------------------------------
    $weightPerCell =~ m|\D| and die "$error: weight-per-cell must be an integer\n";
}

# main execution block
sub svtools_createHMM10X {
    (@globs) = @_;
    $globs[0] and die "$error: createHMM10X does not take any input\n";

    # initialize
    setOptions_createHMM10X();    
    print STDERR "$utility createHMM10X: " . getTime(), "\n";

    # aggregate and normalize data
    print STDERR "$utility createHMM10X: establishing HMM for bin weight $weightPerCell\n";
    $ENV{LIB_DIR}    = $libDir;
    $ENV{HMM_FILE}   = "$libDir/HMM10X/HMM10X.$weightPerCell.RData";
    $ENV{BIN_WEIGHT} = $weightPerCell;
    system("Rscript $libDir/createHMM10X.R") and die "$error: createHMM10X.R returned an error\n";  
}

1;
