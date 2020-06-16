use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            %binFileCols);
my (@globs,
    #$excludeFile,
    $outFile, $copyNum);

# manage options
sub setOptions_exclude {
    #setOptionValue(\$excludeFile,  'exclude-file');
    setOptionValue(\$outFile,      'output-file');
    setOptionValue(\$copyNum,      'copy-number',   5);
    #-------------------------------------
    #$excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");   
}

# main execution block
sub svtools_exclude {
    (@globs) = @_;
    $globs[0] or die "$error: no input file specified\n";
    $globs[1] and die "$error: only one input file allowed for command 'exclude'; ";
    $globs[0] =~ m/svtools\.call\.bins\..+\.bgz/ or
        die "$error: $globs[0] is not a valid call.bins file name\n";

    # initialize
    setOptions_exclude();    
    print STDERR "$utility exclude: " . getTime(), "\n";
    openOutputStream($outFile, \$outH);
    openInputStream($globs[0], \$inH, $FALSE, $FALSE, $TRUE);
   
    # set reference statistics
    print STDERR "$utility exclude: finding consecutive bins with copy number >= $copyNum\n";
    our $binChromCol  = $binFileCols{chrom};    
    our $binStartCol  = $binFileCols{start};
    our $binEndCol    = $binFileCols{end};
    our $binCNCol     = $binFileCols{cnRaw};   
    our ($prevChrom, $sumLen, @excBins) = ("", 0);
    while (my $line = <$inH>){
        chomp $line;
        my @bin = split("\t", $line);
        if ($prevChrom ne $bin[$binChromCol]) {
            $prevChrom and @excBins and commitExclusion();
            $prevChrom = $bin[$binChromCol];
        }
        if ($bin[$binCNCol] >= $copyNum) {
            push @excBins, \@bin;
        } else {
            @excBins and commitExclusion();
        }
    } 
    @excBins and commitExclusion();
    sub commitExclusion {
        my $start = $excBins[0][$binStartCol];
        my $end   = $excBins[$#excBins][$binEndCol];
        my $length = $end - $start;
        my $cn = int(mean(map {$$_[$binCNCol]} @excBins) + 0.5);
        print $outH join("\t", $prevChrom, $start, $end, "svtools.exclude.$length($cn)", 0, "."), "\n";
        @excBins = ();
        $sumLen += $length;
    }    

    # finish up  
    closeHandles($inH, $outH);
    $sumLen = commify($sumLen);
    print STDERR "$utility exclude: $sumLen bp excluded\n";
}

1;
