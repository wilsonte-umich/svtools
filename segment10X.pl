use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            %prbCol);
our ($refineFile,
     $tmpDir, $maxMem,
     $arrayDir, $arrayFmt, $project, $nameCol, $sample,
     $extZyg, $lrrLim, $badPrbFreq, $transProb, $preservation,
     $outputDir, $plotDir);
our ($inputDir, $valType, %chroms);
my ($mdlH, $datH, %mdl);
my ($nMdl, $nArrIn, $nArrMatch) = (0, 0, 0);

# manage options
sub setOptions_segment {
    setOptionValue(\$refineFile,  'refine-file');
    setOptionValue(\$arrayDir,   'array-dir');
    setOptionValue(\$arrayFmt,   'array-format', 'Illumina');    
    setOptionValue(\$project,    'project');
    setOptionValue(\$nameCol,    'name-column',  'DNA_ID');
    setOptionValue(\$sample,     'sample');
    setOptionValue(\$extZyg,       'extreme-zyg',     0.85);
    setOptionValue(\$lrrLim,       'lrr-lim',         '-1.5,1.0');
    setOptionValue(\$badPrbFreq,   'bad-probe-freq',  1e-3);
    setOptionValue(\$transProb,    'transition-prob', 1e-4);
    setOptionValue(\$preservation, 'preservation',    '0.99,0.95,0.9,0.8');
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',         1000000000);  
    #-------------------------------------
    -e $refineFile or die "$error: file not found: $refineFile\n";
    -d $arrayDir or die "$error: directory not found: $arrayDir\n";   
    -d $outputDir or die "$error: directory not found: $outputDir\n";       
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem =~ m|\D| and die "$error: max-mem must be an integer number of bytes\n";
}
# main execution block
sub msvtools_segment {

    # initialize
    setOptions_segment();    
    print STDERR "$utility segment: " . getTime(), "\n";
    
    # load the trained model CN (and read count) for each probe
    loadProbeModel();
    
    # load probe data for target sample
    $inputDir = $arrayDir;
    my $datFile = getOutFile('data', $sample, 'gz');
    openOutputStream($datFile, \$datH, $TRUE);
    print $datH join("\t", qw(
        CHROM START POS PRB_NAME IS_NA STRAND
        RC LRR BAF 
        AS_IN CN_IN INF_IN
        LRR_MED BAF_MED
    )), "\n";
    indexIllumina(\&getProbeData); # not sorted yet, occurs in R at HMM
    closeHandles($datH);
    print STDERR "    $nMdl probes in input model\n";
    print STDERR "    $nArrIn probes in array\n";
    print STDERR "    $nArrMatch probes in array matched the model\n";
    
    # execute the HMM via R
    %mdl = ();    
    segmentArray($sample, $datFile, "$outputDir/plots", "__NULL__",
                 $extZyg, $lrrLim, $badPrbFreq, $transProb, $preservation);
}

sub loadProbeModel {
    print STDERR "$utility segment: loading model: $refineFile\n";
    openInputStream($refineFile, \$mdlH, undef, undef, $TRUE);
    my $header = <$mdlH>;
    while (my $prb = <$mdlH>){
        $nMdl++;        
        chomp $prb;
        my @f = split("\t", $prb);
        $mdl{$f[$prbCol{PRB_NAME}]} = [map { $f[$prbCol{$_}] } qw(AS_OUT CN_OUT INF_OUT
                                                                  LRR BAF)];
        $chroms{$f[$prbCol{CHROM}]}++;
    }
    closeHandles($mdlH);   
}

sub getProbeData {
    my ($f, $hdr) = @_;
    $nArrIn++;
    my $prb = $$f[$$hdr{PRB_NAME}];
    $mdl{$prb} or return;
    $nArrMatch++;
    print $datH join("\t",
        $$f[$$hdr{CHROM}], $$f[$$hdr{POS}]-1, $$f[$$hdr{POS}], $prb, -1, '+',
        0, $$f[$$hdr{LRR}], $$f[$$hdr{BAF}], # this sample's values
        @{$mdl{$prb}}
    ), "\n";
}

1;
