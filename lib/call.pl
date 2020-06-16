use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            %binCols %binColOffsets $nBinLeadCol $nCnvLeadCol);
my (@globs, $group, $delete, $maxMem,
    $reference, $maxCN, $chroms, $minBins);
my (%chroms, @chroms,
    $refChrom, $refStart, $refEnd, $refPloidy,
    $refMedian, $refMean, $refStdev);
my ($gZip, $stepsAfter, $stepsBefore);
my (@samples, $nSamples, @weights);
my (@bins, @cnvs);
my ($binWidth, $minBinSize, $maxBinSize, $maxBinI);
my ($nSmpCNVs, $nGrpCNVs) = (0) x 100;
my $cnvId = 0;

# manage options
sub setOptions_call {
    setOptionValue(\$group,     'group');
    setOptionValue(\$inputDir,  'input-dir');
    setOptionValue(\$delete,    'delete');
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$tmpDir,    'tmp-dir',         '/tmp');
    setOptionValue(\$maxMem,    'max-mem',         1000000000);
    setOptionValue(\$reference, 'reference');
    setOptionValue(\$maxCN,     'max-copy-number', 4);
    setOptionValue(\$chroms,    'chroms');
    setOptionValue(\$minBins,   'min-bins',        4); 
    #-------------------------------------
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    $maxCN  =~ m|\D| and die "$error: option --max-copy-number must be an integer number\n";
    $minBins  =~ m|\D| and die "$error: option --min-bins must be an integer number of bins\n";
    $reference =~ m|(\w+):(\d+)-(\d+)/(\d+)| or die "malformed reference, expected format 'chrom:start-end/copy_number'\n";
    ($refChrom, $refStart, $refEnd, $refPloidy) = ($1, $2, $3, $4);
    %chroms = map { $_ => 1} split(",", $chroms);
    @chroms = sort {$a cmp $b} keys %chroms;
}

# main execution block
sub svtools_call {
    # take no files as input, collected based on group and inputDir

    # initialize
    setOptions_call();    
    print STDERR "$utility call: " . getTime(), "\n";
    my $smpFile = getOutFile("samples", $group, "txt");
    my $binFile = getOutFile("bins",    $group, "bgz");
    my $cnvFile = getOutFile("cnvs",    $group, "bgz");
    openOutputStream($smpFile, \my $smpH);
    openOutputStream($binFile, \my $binH, undef, undef,
                     "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    openOutputStream($cnvFile, \my $cnvH, undef, undef,
                     "Rscript $libDir/call.R $nCnvLeadCol | sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");

    # set reference statistics
    print STDERR "$utility call: determining statistics for reference region\n";
    setRefStats();
    initializeCopyNumberHMM();

    # solve for bin copy number states by HMM
    # thread data for runs of unusual copy number, i.e. CNVs
    print STDERR "$utility call: solving for copy number states\n";
    foreach my $binFile(glob(getInFile("bin", "bins", "$group.*", "gz"))){
        my $chrom = loadBins($binFile);
        $chroms{$chrom} or next;
        print STDERR "  $chrom\n";
        @cnvs = ();
        solveCopyNumberHMM();
        unless(@samples == 1){
            foreach my $sampleI(0..$#samples){ solveCnvHMM($sampleI) }        
            foreach my $sampleI(0..$#samples){ findSampleCNVs($sampleI) }
            findGroupCNVs();            
            foreach my $cnv(@cnvs){
                print $cnvH join("\t", $chrom, @$cnv), "\n";
            }            
        }
        foreach my $bin(@bins){
            print $binH join("\t", $chrom, @$bin), "\n";
        }
        $delete and unlink $binFile;
    }
    
    # finish up
    unlink $ENV{eProbFile_CN};
    print STDERR "$utility call: finalizing results\n";
    print $smpH join("\t", @samples), "\n";    
    closeHandles($smpH, $binH, $cnvH);
    system("tabix -p bed $binFile");
    system("tabix -p bed $cnvFile");
    print STDERR "$utility call: $nSmpCNVs CNVs called in individual samples\n";
    print STDERR "$utility call: $nGrpCNVs CNV regions called for group $group\n";
}

# load the bins for a specific chromosome
sub loadBins {
    my ($binFile) = @_;
    $binFile =~ m/svtools\.bin\.bins\.$group\.(\w+)/ or
        die "$utility bin error: $binFile is not a valid bins file name for $group\n";    
    my $chrom = $1;
    openInputStream($binFile, \$inH, $FALSE, $FALSE, $TRUE);
    my $line = <$inH>;
    chomp $line;
    $line =~ m/^#(.+)/ or die "$utility call error: no header line in bin file: $binFile\n";    
    my %tags = getFileHeader($1, 'bin', qw(group samples weights)); 
    @samples = split(",", $tags{samples});
    $nSamples = @samples;
    @weights = split(",", $tags{weights});
    @bins = ();    
    while (my $line = <$inH>){
        chomp $line;
        push @bins, [split("\t", $line)];
    } 
    closeHandles($inH);
    return $chrom;
}

# get information on the main peak of adjusted bin sizes (for CN assessment)
sub setRefStats {
    loadBins(getInFile("bin", "bins", "$group.$refChrom", "gz"));
    my $startCol = $binCols{start};
    my $endCol   = $binCols{end};
    my $lenCol   = $binCols{adjLen};
    my @sizes;
    # restrict attention to the instructed reference region
    foreach my $bin(@bins){ 
        $$bin[$startCol] >= $refStart or next;
        $$bin[$endCol]   <= $refEnd   or next;
        push @sizes, $$bin[$lenCol];
    }
    # find the sweet spot of the bin size distribution
    my $median = median(@sizes);
    my $min = $median - $median / 2;
    my $max = $median + $median / 2;
    # calculate stats on the bin size peak, assuming Normal (general valid)
    my @peak;    
    foreach my $size(@sizes){
            $size >= $min
        and $size <= $max
        and push @peak, $size;
    }
    ($refMedian, $refMean, $refStdev) = (median(@peak), stdev(@peak));
    my $rmean  = int($refMean + 0.5);
    my $rstdev = int($refStdev + 0.5);
    print STDERR "$utility call: $rmean +/- $rstdev reference bin sizes => CN$refPloidy\n"
}

# initialize the copy number HMM for genome based on reference region
sub initializeCopyNumberHMM {
    # set the boundaries for the bin-width HMM
    my $nSD   = 3.5;  # helps determine HMM boundary limits
    my $nBins = 100; # number of divisions, i.e. observation states, in the HMM (actually get one more)
    $maxBinI = $nBins + 2;
    my ($mean1, $stdev1) = adjustSizeStats(1);       # CN=1, largest varBin size
    my ($meanX, $stdevX) = adjustSizeStats($maxCN); # smallest varBin size
    my $maxSize = $mean1 + $nSD * $stdev1; # largest and smallest varBin sizes to allow in HMM
    my $minSize = $meanX - $nSD * $stdevX;
    $binWidth = int(($maxSize - $minSize) / $nBins + 0.5); # width of each bin in unit varBin size
    $maxBinSize = getBinSize($maxSize); # rounded bin values at the limit bins in unit varBin length
    $minBinSize = getBinSize($minSize);
    # generate a set of copy number emission probabilities for the genome based on reference
    $ENV{eProbFile_CN} = "$tmpDir/svtools.$group.emiss.prob.".int(rand(1e6)).".txt";
    $ENV{refMean}    = $refMean;
    $ENV{refStdev}   = $refStdev;
    $ENV{refPloidy}  = $refPloidy;
    $ENV{binWidth}   = $binWidth;
    $ENV{maxBinSize} = $maxBinSize;
    $ENV{minBinSize} = $minBinSize;
    $ENV{maxCN}      = $maxCN;
    system("Rscript $libDir/set_CN_emissions.R") == 0
        or die "error setting CN emission probabilities\n";     
}

# solve for the modal copy number across all samples
sub solveCopyNumberHMM { 
    my $ajdLenCol = $binCols{adjLen};
    my $cnRawCol  = $binCols{cnRaw};
    my $cnHMMCol  = $binCols{cnHMM};

    # convert bin data to observation indices, correlated to bin indices
    my $tmpFile = "$tmpDir/svtools.$group.CN.segment.".int(rand(1e6)).".data";  
    open(my $tmpH, ">", $tmpFile) or die "could not open $tmpFile for writing: $!\n";
    foreach my $binI(0..$#bins){ 
        print $tmpH join("\t", getBinI($bins[$binI][$ajdLenCol]), $binI), "\n"; 
    }
    close $tmpH;
    
    # convert resulting CN state indices to copy number
    open(my $hmmH, "-|", "cat $tmpFile | segment -e $ENV{eProbFile_CN} -z 0.1 -p 0.95")
        or die "could not open CN segment stream: $!\n";
    while (my $line = <$hmmH>) {
        chomp $line;
        my ($CN, $binI) = split("\t", $line);
        $CN > $maxCN and $CN = $maxCN;
        my $bin = $bins[$binI];
        $$bin[$cnRawCol] = roundCount($refMean / $$bin[$ajdLenCol] * $refPloidy, 100);         
        $$bin[$cnHMMCol] = $CN; 
    }
    close $hmmH;
    unlink $tmpFile;    
}

# solve for the copy number change of each sample
sub solveCnvHMM { 
    my ($sampleI) = @_;
    my $sample     = $samples[$sampleI];
    $ENV{weight}   = $weights[$sampleI]; # sample level data sent to R      
    my $medCovCol  = $binCols{medCov};    
    my $cnRawCol   = $binCols{cnRaw};
    my $cnHMMCol   = $binCols{cnHMM};
    my $smpCovRawCol = $nBinLeadCol + $nSamples*$binColOffsets{smpCovRaw} + $sampleI;
    my $smpCovNrmCol = $nBinLeadCol + $nSamples*$binColOffsets{smpCovNrm} + $sampleI;
    my $smpDCNRawCol = $nBinLeadCol + $nSamples*$binColOffsets{smpDCNRaw} + $sampleI;
    my $smpDCNHMMCol = $nBinLeadCol + $nSamples*$binColOffsets{smpDCNHMM} + $sampleI;

    # set the tmp files
    $ENV{dataFile}      = "$tmpDir/svtools.$group.CNV.segment.".int(rand(1e6)).".data";  
    $ENV{eProbFile_CNV} = "$tmpDir/svtools.$group.CNV.emiss.prob.".int(rand(1e6)).".txt";    

    # extract bin-level data for R, correlated to bin indices
    open(my $tmpH, ">", $ENV{dataFile}) or die "could not open $ENV{dataFile} for writing: $!\n";
    foreach my $binI(0..$#bins){
        print $tmpH join("\t", (map {$bins[$binI][$_]}($medCovCol, $cnHMMCol, $smpCovRawCol)), $binI), "\n";
    }
    close $tmpH;

    # calculate the bin-specific emission probabilities
    system("Rscript $libDir/CNV_HMM.R") == 0
        or die "error setting CNV emission probabilities\n";

    # convert resulting cnDelta state indices to copy number changes
    open(my $hmmH, "-|", "cat $ENV{eProbFile_CNV} | segment -z 0.95 -p 0.95 -w")
        or die "could not open CNV segment stream: $!\n";
    while (my $line = <$hmmH>) {
        chomp $line;
        my ($dCN, $binI) = split("\t", $line);
        my $bin = $bins[$binI];
        $$bin[$smpDCNRawCol] =
            roundCount(($$bin[$smpCovNrmCol] - $$bin[$medCovCol]) / $$bin[$medCovCol] * $$bin[$cnRawCol], 100);
        $$bin[$smpDCNHMMCol] = $dCN - 2;
    } 
    close $hmmH;
    unlink $ENV{dataFile};
    unlink $ENV{eProbFile_CNV};
}

# find runs of CNV bins per sample
sub findSampleCNVs {
    our ($sampleI) = @_;
    our $binStartCol  = $binCols{start};    
    our $binEndCol    = $binCols{end};
    our $cnvCallsCol  = $binCols{cnvCalls};
    our $smpDCNHMMCol = $nBinLeadCol + $nSamples*$binColOffsets{smpDCNHMM}  + $sampleI;
    our $smpCNVCallCol= $nBinLeadCol + $nSamples*$binColOffsets{smpCNVCall} + $sampleI;
    our ($inRun, @cnvBins);
    foreach my $bin(@bins){
        $$bin[$cnvCallsCol] = $$bin[$cnvCallsCol] eq 'na' ? 0 : $$bin[$cnvCallsCol];
        $$bin[$smpCNVCallCol] = 0;
        my $cnv = $$bin[$smpDCNHMMCol] < 0 ? 'loss' : ($$bin[$smpDCNHMMCol] > 0 ? 'gain' : '');   
        defined $inRun and $cnv ne $inRun and commitSampleCNV();  
        if ($cnv) {
            $inRun = $cnv;
            push @cnvBins, $bin;
        }
    }
    $inRun and commitSampleCNV();
    sub commitSampleCNV {
        my $nCnvBins = @cnvBins; 
        if($nCnvBins >= $minBins) {
            $nSmpCNVs++;
            my $start = $cnvBins[0][$binStartCol];
            my $end   = $cnvBins[$#cnvBins][$binEndCol];
            $cnvId++;
            push @cnvs, [$start, $end, $cnvId, $nCnvBins, '.',
                         $inRun, $end-$start, $samples[$sampleI], $sampleI,
                         describeCNVBins(\@cnvBins), describeCNVSamples(\@cnvBins)];
            my $cnvCall = $inRun eq 'loss' ? -1 : 1; # signed boolean to represent gain vs. loss in this bin for this sample
            foreach my $cnvBin(@cnvBins){
                $$cnvBin[$cnvCallsCol]++;
                $$cnvBin[$smpCNVCallCol] = $cnvCall;
            }
        }
        $inRun = undef;
        @cnvBins = ();
    }   
}

  #chrom start end cnvId nBins strand
  #cnvType cnvSize sample sampleIndex
  #uniqueness maxNSamples binsList

# find runs of CNV bins for the group, i.e. collapse the per-sample CNV calls
sub findGroupCNVs {
    our $binStartCol  = $binCols{start};    
    our $binEndCol    = $binCols{end};
    our $cnvCallsCol  = $binCols{cnvCalls};
    our @cnvBins;
    foreach my $bin(@bins){
        @cnvBins and !$$bin[$cnvCallsCol] and commitGroupCNV();     
        $$bin[$cnvCallsCol] and push @cnvBins, $bin;
    }
    @cnvBins and commitGroupCNV();
    sub commitGroupCNV {
        $nGrpCNVs++;
        my $nCnvBins = @cnvBins; 
        my $start = $cnvBins[0][$binStartCol];
        my $end   = $cnvBins[$#cnvBins][$binEndCol];
        $cnvId++;
        push @cnvs, [$start, $end, $cnvId, $nCnvBins, '.',
                     'cnv', $end-$start, 'group', -1,
                     describeCNVBins(\@cnvBins), describeCNVSamples(\@cnvBins)];
        @cnvBins = ();
    }  
}

# calculate aggregate metrics to determine the uniqueness of a run of CNV bins
sub describeCNVBins {
    my ($cnvBins) = @_; 
    my @list;
    foreach my $cnvBin(@$cnvBins){
        my $nSmp = 0; 
        foreach my $sampleI(0..$#samples){
            my $smpDCNHMMCol = $nBinLeadCol + $nSamples*$binColOffsets{smpDCNHMM} + $sampleI;
            $nSmp += ($$cnvBin[$smpDCNHMMCol] ?  1 : 0); # boolean that HMM called CN deviation in bin for this sample,
        }                                                # based on each bin's HMM, NOT considering minBins calling threshold
        push @list, $nSmp; # ordered collection of nSamples with CNV in each bin
    }
    my $uniqueness = roundCount(mean(@list), 100);
    my ($maxNSmp, $prevNSmp, $nBins, $list) = (0, 0, 0, '');
    foreach my $nSmp(@list){
        $maxNSmp >= $nSmp or $maxNSmp = $nSmp;
        if($prevNSmp and $prevNSmp != $nSmp) {
            $list .= "$prevNSmp($nBins)";
            $nBins = 0;
        }
        $prevNSmp = $nSmp;
        $nBins++;
    }
    $list .= "$prevNSmp($nBins)";
    return ($uniqueness, $maxNSmp, $list);
}

# aggregate the delta copy number over all bins in a CNV for all samples
sub describeCNVSamples {
    my ($cnvBins) = @_;
    my (@means, @dCNs);
    foreach my $sampleI(0..$#samples){
        my $smpDCNRawCol = $nBinLeadCol + $nSamples*$binColOffsets{smpDCNRaw} + $sampleI;
        my @smpDCNRaw = map {$$_[$smpDCNRawCol]} @$cnvBins;
        push @means, roundCount(mean(@smpDCNRaw), 100);
        push @dCNs, join(",", @smpDCNRaw);
    }
    return (@means, @dCNs);  
}

# bin size subs
sub getBinSize {
    my ($varBinSize) = @_;
    return int($varBinSize / $binWidth + 0.5) * $binWidth; 
}
sub getBinI {
    my ($varBinSize) = @_;
    my $sizeBin = getBinSize($varBinSize);
    if ($sizeBin < $minBinSize ) {
        return 0;
    } elsif($sizeBin > $maxBinSize){
        return $maxBinI;
    } else {
        return ($sizeBin - $minBinSize) / $binWidth + 1;
    }
}
sub adjustSizeStats {
    my ($CN) = @_;
    return ($refMean  * $refPloidy / $CN, 
            $refStdev * $refPloidy / $CN); 
}

1;


## WHAT IS THIS FOR???
#sub setBinStats {
#    my $lenCol = $binCols{adjLen};    
#    my $nSD = 2.5;
#    foreach my $chrom(@chroms){
#        $skipChrom{$chrom} and next;
#        my %sizes;    
#        foreach my $CN(1..$refPloidy+1){
#            my ($cnMean, $cnStdev) = adjustSizeStats($CN);
#            my $min = $refMean - $refStdev * $nSD;
#            my $max = $refMean + $refStdev * $nSD;
#            foreach my $bin(@{$bins{$chrom}}){
#                    $$bin[$lenCol] >= $min
#                and $$bin[$lenCol] <= $max
#                and push @{$sizes{$CN}}, $$bin[$lenCol];
#            }              
#        }
#        my @nullChrom = (0, 0, 0);
#        my @availCNs = keys %sizes;
#        @availCNs or return @nullChrom; 
#        my @CNs = sort { @{$sizes{$b}} <=> @{$sizes{$a}} } @availCNs;    
#        my $CN = $CNs[0];
#        @{$sizes{$CN}} > 100 or return @nullChrom; 
#        my @peak = @{$sizes{$CN}};
#        $binStats{$chrom} = [median(@peak), stdev(@peak)];# median, mean, stdev
#    }
#}