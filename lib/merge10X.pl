use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH
            $inputDir $tmpDir $outputDir $plotDir);
require "$libDir/exclude_subs.pl";
my (@globs, $sample, $excludeFile);
my ($nCells, $nBins, $nChrom) = (0) x 10;

# manage options
sub setOptions_merge10X {
    setOptionValue(\$inputDir,  'input-dir');
    setOptionValue(\$sample,    'sample');
    setOptionValue(\$excludeFile,    'exclude-file', '');
    setOptionValue(\$outputDir, 'output-dir');
    #-------------------------------------
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    $outputDir or $outputDir = $inputDir;
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
    $plotDir = "$outputDir/plots";
    mkdir $plotDir;
    -d $plotDir or die "$error: failed to create $plotDir\n";
}

# main execution block
sub svtools_merge10X {
    (@globs) = @_;
    $globs[0] and die "$error: merge10X does not take any input (bin files are read from disk)\n";

    # initialize
    setOptions_merge10X();    
    print STDERR "$utility merge10X: " . getTime(), "\n";
    initializeExclude($excludeFile);
    
    # extract common header
    print STDERR "$utility merge10X: extracting header\n";  
    my $inFiles = getInFile('bin10X', 'bins', "$sample.*", 'gz');
    my @inFiles = glob($inFiles);
    openInputStream($inFiles[0], \my $hdrH, $TRUE);
    my $header = <$hdrH>;
    chomp $header;
    $header =~ m/^#/ or $header = "#".$header;
    my ($chrom, $start, $end, $rest) = split("\t", $header, 4);
    $header = join("\t", $chrom, $start, $end, "nGapBases", $rest);
    close $hdrH;
    
    # merge chromosomes and count
    print STDERR "$utility merge10X: merging bins over all chromosomes\n";
    my $outFile = getOutFile("bins", $sample, "bgz");
    my $stepsBefore = "sort -S 1G -k1,1 -k2,2n -k3,3n | ".
                      "awk 'BEGIN{print \"$header\"}{print \$0}' | ".
                      "bgzip -c";
    openOutputStream($outFile, \$outH, $FALSE, "50M", $stepsBefore);
    foreach my $file(@inFiles){
        $nChrom++;
        openInputStream($file, \$inH, $TRUE);
        my $header = <$inH>;
        my @header = split("\t", $header);
        my $fileNCells = @header - 4;
        if($nCells){
            $nCells == $fileNCells or die "$error: input files have inconsistent numbers of cells\n";               
        } else {
            $nCells = $fileNCells;
        }
        while (my $line = <$inH>){
            $nBins++;              
            my ($chrom, $start, $end, $rest) = split("\t", $line, 4);
            print $outH join("\t", $chrom, $start, $end,
                                   countExcluded($chrom, $start, $end), $rest);
        }
        closeHandles($inH);
    }
    closeHandles($outH);
    system("tabix -p bed $outFile") and die "$error: tabix returned an error\n";
    
    # TODO: handle deletion of chromosome bin files    
    
    # report counts
    print STDERR "\n";
    reportCount($nCells,          "true cells");
    reportCount($nChrom,          "chromosomes");
    reportCount($nBins,           "bins over all chromosomes");    
    print STDERR "\n";

    # aggregate and normalize data
    print STDERR "$utility merge10X: normalizing data and plotting summaries\n";
    $ENV{LIB_DIR}    = $libDir;
    $ENV{IN_DIR}     = $inputDir;
    $ENV{OUT_DIR}    = $outputDir;
    $ENV{PLOT_DIR}   = $plotDir;
    $ENV{SAMPLE}     = $sample;
    system("Rscript $libDir/merge10X.R") and die "$error: merge10X.R returned an error\n";  
}

1;
