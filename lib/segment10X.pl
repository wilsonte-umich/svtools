use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH
            $inputDir $tmpDir $outputDir $plotDir);
my (@globs, $sample, $chrom,
    $transProb, $preservation, $weightPerCell);
my ($nCells, $nBins, $nChrom) = (0) x 10;

# manage options
sub setOptions_segment10X {
    setOptionValue(\$inputDir,  'input-dir');
    setOptionValue(\$sample,    'sample');
    setOptionValue(\$chrom,     'chrom');
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$transProb, 'transition-prob', 1e-4);
    setOptionValue(\$preservation, 'preservation', "0.99,0.95,0.9,0.8");
    setOptionValue(\$weightPerCell, 'weight-per-cell', 10);     
    #-------------------------------------
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    $outputDir or $outputDir = $inputDir;
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
    $plotDir = "$outputDir/plots";
    mkdir $plotDir;
    -d $plotDir or die "$error: failed to create $plotDir\n";
    $weightPerCell =~ m|\D| and die "$error: weight-per-cell must be an integer\n";
}

# main execution block
sub svtools_segment10X {
    (@globs) = @_;
    $globs[0] and die "$error: segment10X does not take any input (bin files are read from disk)\n";

    # initialize
    setOptions_segment10X();    
    print STDERR "$utility segment10X: " . getTime(), "\n";

    # aggregate and normalize data
    print STDERR "$utility segment10X: segmenting chromosome $chrom\n";
    $ENV{LIB_DIR}    = $libDir;
    $ENV{HMM_FILE}   = "$libDir/HMM10X/HMM10X.$weightPerCell.RData";
    $ENV{BIN_WEIGHT} = $weightPerCell;
    $ENV{BINS_FILE}  = getInFile('merge10X', "bins", $sample, "bgz");
    $ENV{CNC_FILE}   = getInFile('merge10X', "CNC",  $sample, "bgz");
    $ENV{TRANS_PROB} = $transProb;
    $ENV{PRESERVATION} = $preservation;
    $ENV{PLOT_DIR}   = $plotDir;
    $ENV{SAMPLE}     = $sample;
    $ENV{CHROM}      = $chrom;
    $ENV{CALL1_FILE} = getOutFile("call_single", "$sample.$chrom", "bgz");
    $ENV{PLOT_DIR} = $plotDir;

    
    system("Rscript $libDir/segment10X.R") and die "$error: segment10X.R returned an error\n";  
}

1;
