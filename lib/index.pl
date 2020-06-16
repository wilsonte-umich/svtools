use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $_gap $_split $_clip
            %idxSubs $idxH $libKey);
require "$libDir/index_common.pl";
our ($TRUE, $FALSE) = (1, undef);
my (@globs, $group, $maxMem, $minTLen, $useMask);
our ($tmpDir, $outputDir, %mask);
my ($inH);

# manage options
sub setOptions_index {  
    setOptionValue(\$group,     'group');
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$maxMem,    'max-mem',       1000000000);    
    setOptionValue(\$minTLen,   'min-tLen',      '1.5/-1/0');
    setOptionValue(\$useMask,   'use-mask');
    #-------------------------------------
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
}

# main execution block
sub svtools_index {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";
    
    # initialize
    setOptions_index();    
    print STDERR "$utility index: " . getTime(), "\n";

    # open output stream
    my $outFile = getOutFile('anomalies', $group, "bgz");
    openOutputStream($outFile, \$idxH, undef, undef,
                     "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    
    # do the work
    processInputFiles(\$inH, \@globs, $TRUE, $FALSE, $TRUE, \&indexExtractFile);
    closeHandles($idxH);
    system("tabix -p bed $outFile");
}

#========================================================================
# index 'extract'
#------------------------------------------------------------------------
sub indexExtractFile {
    my ($file) = @_;
    # read the extract header
    my %minTLen;
    my $line = <$inH>;
    chomp $line;
    $line =~ m/^#(.+)/ or die "$error: missing header in $file\n";
    my %tags = getFileHeader($1, 'extract', qw(sample library maxTLen medTLen));
    $libKey  = "$tags{sample}.$tags{library}";
    setMinTLen('index', $libKey, $minTLen, \%minTLen, \%tags);    
    # process the anomalies
    while (my $line = <$inH>) {
        chomp $line;
        my ($pairType, $rest) = split("\t", $line, 2);
        checkTLen($pairType, $rest, \%minTLen) or next;
        &{$idxSubs{$pairType}}(split("\t", $rest)); # TODO: stop including clips?
    }
    closeHandles($inH);
}
#========================================================================

1;
