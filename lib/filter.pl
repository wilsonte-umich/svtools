use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir $slurp
            $TRUE $FALSE
            $tmpH $inH $outH
            $inputDir $tmpDir $outputDir $plotDir
            %gntCols %poolCols);
my (@globs, $group, $maxMem,
    $minReads, $maxReads, $minFracAlt, $maxFracAlt);
my %c = %poolCols;

# manage options
sub setOptions_filter {
    setOptionValue(\$group,      'group');
    setOptionValue(\$inputDir,   'input-dir');
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$minReads,   'min-reads', 10);
    setOptionValue(\$maxReads,   'max-reads', 100);
    setOptionValue(\$minFracAlt, 'min-frac-alt', 0.3);
    setOptionValue(\$maxFracAlt, 'max-frac-alt', 0.7);
    #-------------------------------------
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    $minReads =~ m|\D| and die "$error: min-reads must be an integer\n";
    $maxReads =~ m|\D| and die "$error: max-reads must be an integer\n";    
}

# main execution block
sub svtools_filter {

    # initialize
    setOptions_filter();    
    print STDERR "$utility filter: " . getTime(), "\n";
    
    # filter heterozygous SNPs
    my $poolFile = getInFile('pool', 'snps', $group, "bgz");
    openInputStream("gunzip -c $poolFile", \$inH);    
	my $fltFile = getOutFile("snps", $group, "bgz");
    openOutputStream($fltFile, \$outH, undef, undef,
            "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    my $nKept = 0;
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        if ($f[$c{nReads}] >= $minReads and
            $f[$c{nReads}] <= $maxReads and
            $f[$c{nAlt}]/$f[$c{nReads}] >= $minFracAlt and
            $f[$c{nAlt}]/$f[$c{nReads}] <= $maxFracAlt ) {
            print $outH "$line\n";
            $nKept++;
        }  
    }
	closeHandles($inH, $outH);    

    # finish up
    print STDERR commify($nKept)," SNPs kept as heterozygous\n";
    system("tabix -p bed $fltFile");
}

1;
