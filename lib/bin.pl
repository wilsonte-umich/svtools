use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir $slurp
            $TRUE $FALSE
            $inH $tmpH $outH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            %binCols %binColOffsets $nBinLeadCol);
my (@globs, $group, $gapFile, $excludeFile, $maxMem,
    $chrom, $chroms, $binCount, $refSamples);
my ($awk, $gZip, $stepsAfter, $stepsBefore);
my (@refSamples, %refSamples, @refFiles);
my (%chroms, @chroms, %chromSizes);
my (%libStats, %smpStats, %files, @samples, $nSamples, @weights);
my (%breaks_XXX, @bins, @counts);
my ($nXXX) = (0) x 100;

# manage options
sub setOptions_bin {
    setOptionValue(\$group,      'group');
    setOptionValue(\$gapFile,    'gap-file');
    setOptionValue(\$excludeFile,'exclude-file');
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$tmpDir,     'tmp-dir',         '/tmp');
    setOptionValue(\$maxMem,     'max-mem',         1000000000);
    setOptionValue(\$chrom,      'chrom');
    setOptionValue(\$chroms,     'chroms');
    setOptionValue(\$binCount,   'bin-count',       100);
    setOptionValue(\$refSamples, 'ref-samples');
    #-------------------------------------
    $gapFile and (-e $gapFile or die "$error: file not found: $gapFile\n");
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";   
    $binCount  =~ m|\D| and die "$error: option --bin-count must be an integer number of pairs\n";
    $awk = "awk '\$1==\"$chrom\"'";
    %chroms = map { $_ => 1} ($chrom, split(",", $chroms));
    @chroms = sort {$a cmp $b} keys %chroms;
    $refSamples and @refSamples = split(",", $refSamples);
    %refSamples = map { $_ => 1} @refSamples;
}

# main execution block
sub svtools_bin {
    (@globs) = @_;
    $globs[0] or die "$error: no input file or stream specified\n";
    
    # initialize
    setOptions_bin();    
    print STDERR "$utility bin: " . getTime(), "\n";
    processInputFiles(\$inH, \@globs, $FALSE, '', $TRUE, \&get_map_header);
    setWeights();    
    
    # define variable width bins based on all samples (or all reference samples)
    print STDERR "$utility bin: setting ends of variable width bins on $chrom\n";
    my $refFiles = $refSamples ? \@refFiles : \@globs;
    my $tmpFile = getTmpFile('bin', 'break_map_all', 'gz');    
    openOutputStream($tmpFile, \$tmpH, $TRUE, undef, "sort -T $tmpDir -S $maxMem"."b -k1,1n | groupBy -g 1 -c 2 -o sum");    
    processInputFiles(\$inH, $refFiles, $FALSE, $awk, $TRUE, \&loadBreakMap);    
    closeHandles($tmpH);
    openInputStream($tmpFile, \$tmpH);
    setBinEnds();
    closeHandles($tmpH);
    unlink $tmpFile; 
    #%breaks = ();
    
    # count fragments in each bin for each sample
    print STDERR "$utility bin: counting bin coverage per sample\n";
    foreach my $sample(@samples){
        print STDERR "  $sample\n";
        @counts = ();
        processInputFiles(\$inH, $files{$sample}, $FALSE, $awk, $TRUE, \&countSample);        
        map { push @{$bins[$_]}, $counts[$_] || 0 } 0..$#bins;                  
    }
    @counts = ();
    
    # establish aggregate information on bins
    print STDERR "$utility bin: calculating aggegrate bin statistics\n";
    adjustBinSizes();

    # commit the  results
    print STDERR "$utility bin: storing results on disk\n";
    my $binFile = getOutFile("bins", "$group.$chrom", "gz");
    openOutputStream($binFile, \my $binH, $TRUE);
    print $binH join("\t", "#application:svtools",
                            "version:$version",
                            "file_type:bin",
                            "group:$group",
                            "bin-count:$binCount",
                            "chroms:".join(",", @chroms),
                            "samples:".join(",", @samples),
                            "weights:".join(",", @weights)), "\n"; 
    foreach my $bin(@bins){
        print $binH join("\t", @$bin), "\n";
    }
    closeHandles($binH);
    my $nBins = commify(scalar(@bins));
    my $chromSize = commify($chromSizes{$chrom});
    my $medianSize = commify(median(map { $$_[$binCols{adjLen}] } @bins));
    print STDERR "$utility bin: $nBins bins in $chromSize bp of $chrom\n";
    print STDERR "$utility bin: median adjusted bin size = $medianSize bp\n";
}

# parse the map header information about each sample library file
sub get_map_header {
    my ($file) = @_;
    my $line = <$inH>;
    chomp $line;
    $line =~ m/^#(.+)/ or die "$utility bin error: no header line in map file: $file\n";    
    my %tags = getFileHeader($1, 'map', qw(sample library chromSizes chromStats medTLen));
    print STDERR "$utility bin: library = $tags{sample}.$tags{library}\n";
    push @{$files{$tags{sample}}}, $file;
    $refSamples{$tags{sample}} and push @refFiles, $file;    
    foreach my $cs(split(";", $tags{chromSizes})){
        my ($chrom, $size) = split(",", $cs);
        $chromSizes{$chrom} = $size;
    }
    foreach my $cs(split(";", $tags{chromStats})){
        my ($chrom, $nFrags, $avgLen, $avgCov) = split(",", $cs);
        $chroms{$chrom} or next;
        $libStats{$file}{$chrom} = [$nFrags, $avgLen, $avgCov, $tags{medTLen}];
        $smpStats{$tags{sample}}{$chrom} += $nFrags;
        $smpStats{$tags{sample}}{all} += $nFrags;  
    }
    closeHandles($inH);
}

# parse integrated information from the sample headers
sub setWeights {
    @samples = sort {$a cmp $b} keys %smpStats;
    $nSamples = @samples;
    my $medianAll = median(map {$smpStats{$_}{all}} @samples);
    @weights = map {$smpStats{$_}{all} / $medianAll} @samples;
}

# merge the coverage maps of all input samples (or all reference samples)
sub loadBreakMap {
    my ($file) = @_;
    my $avgLen = $libStats{$file}{$chrom}[1]; # working one chromosome at a time 
    while (my $line = <$inH>) {
        $line =~ m/^#/ and next;
        chomp $line; # sum the coverage changes at all break points in all files across all samples
        my ($chrom, $break, $inc) = split("\t", $line);
        print $tmpH join("\t", $break, $inc / $avgLen), "\n";
        #$breaks{$break} 
        #    = ($breaks{$break} || 0)
        #    +  $inc / $avgLen;    
    }
    closeHandles($inH); 
}

# thread the merged breaks to find the optimal bin endpoints based on all samples
sub setBinEnds {
    my $multiplier = $refSamples ? scalar(@refSamples) : scalar(@samples);
    my $binCountAll = $binCount * $multiplier;
    my ($prevBreak, $cov, $binSum, $binStart) = (0) x 10;
    my @trailing = ('.', 'na', 'na', 'na', 'na', 'na', 'na');
    #BREAK: foreach my $break(sort {$a <=> $b} keys %breaks){    
    BREAK: while (my $line = <$tmpH>) {
        chomp $line;
        my ($break, $inc) = split("\t", $line);
        my $spanCount = ($break - $prevBreak) * $cov;
        $binSum += $spanCount;
        while ($binSum >= $binCountAll and $spanCount) {
            my $binEnd = int(0.5 + $prevBreak - 1 +
                             ($break - $prevBreak) * (1 - ($binSum - $binCountAll) / $spanCount));
            push @bins, [$binStart, $binEnd, '.', $binEnd - $binStart + 1, @trailing];
            $binStart = $binEnd + 1;        
            $binSum = ($break - $binStart) * $cov;
            $prevBreak = $binStart;
            $spanCount = $binSum;
        }
        $prevBreak = $break;
        #$cov += $breaks{$break};
        $cov += $inc;
    }
    # finish the chromosome
    if ($binStart <= $chromSizes{$chrom}) {
        push @bins, [$binStart, $chromSizes{$chrom}, '.', $chromSizes{$chrom} - $binStart + 1, @trailing];
    }
    ($gapFile or $excludeFile) and correctForGaps();
}

# correct bin lengths by subtracting out gap+exclusion lengths
sub correctForGaps {
    # write this chromosome's bins to BED for bedtools
    my $tmpFile = "$tmpDir/$group.$chrom.gap.".int(rand(1e6)).".tmp.bed";
    open(my $outH, ">", $tmpFile) or die "could not open $tmpFile for writing: $!\n";
    my @cols = ($binCols{start}, $binCols{end});
    foreach my $bin(@bins){
        print $outH join("\t", $chrom, @$bin[@cols]), "\n";
    }
    close $outH;
    # run bedtools against the gap/exclusion files, counting total gap+exclusion length
    $gapFile or $gapFile = "";
    $excludeFile or $excludeFile = "";
    open(my $inH, "-|",
        "cat $gapFile $excludeFile | cut -f 1-3 | ".
        "bedtools intersect -wao -a $tmpFile -b - | ".
        "awk 'BEGIN{OFS=\"\t\"}{print \$1,\$2,\$3,\$NF}' | ".
        "groupBy -g 1,2,3 -c 4 -o sum | cut -f 4"
    ) or die "could not open correctForGaps stream: $!\n";
    # adjust the actual gap length or each bin by subtracting its gap+exclusion length
    my $i = 0;
    my $lenCol      = $binCols{length};
    my $noGapLenCol = $binCols{noGapLen};
    while (my $gapLen = <$inH>) {
        chomp $gapLen;
        $gapLen or $gapLen = 0;
        $bins[$i][$noGapLenCol] = $bins[$i][$lenCol] ? $bins[$i][$lenCol] - $gapLen : 0;
        $i++;
    }
    close $inH;
    unlink $tmpFile;
}

# thread each sample to find its counts in each bin
sub countSample {
    my ($file) = @_;
    our $binStartCol = $binCols{start};    
    our $binEndCol   = $binCols{end};
    my $avgLen = $libStats{$file}{$chrom}[1];
    our ($prevBreak, $cov, $binSum, $binI) = 
        (0,          0,    0,       0);
    while (my $line = <$inH>) {
        $line =~ m/^#/ and next;
        chomp $line;
        my ($chrom, $break, $inc) = split("\t", $line);
        processBreak($break);
        $prevBreak = $break;
        $cov += $inc / $avgLen;
        $binI > $#bins and last;
    }
    my $lastBreak = $chromSizes{$chrom} + 1;
    $lastBreak > $prevBreak and processBreak($lastBreak);  
    closeHandles($inH);
    sub processBreak { # fill all necessary bins with the distributed count from this break
        my ($break) = @_;    
        my $spanCount = ($break - $prevBreak) * $cov;
        $binSum += $spanCount;
        while ($bins[$binI] and $break > $bins[$binI][$binEndCol] and $break > $prevBreak) {
            my $excess = $spanCount * ($break - $bins[$binI][$binEndCol] - 1) / ($break - $prevBreak);
            my $count = roundCount($binSum - $excess);
            $counts[$binI] += $count; # keep adding each file for a sample to its cumulative counts
            $binSum = $excess;
            $binI++;
            $binI > $#bins and last;
            $prevBreak = $bins[$binI][$binStartCol];
            $spanCount = $binSum;   
        }   
    }   
}

# calculated the adjusted bin sizes, i.e. correcting for the effect of CNV-deviant samples
sub adjustBinSizes{
    my $maxSampleI  = $#samples;
    my $medCovCol   = $binCols{medCov};
    my $adjLenCol   = $binCols{adjLen};
    my $noGapLenCol = $binCols{noGapLen};
    foreach my $bin(@bins){
        defined $$bin[$nBinLeadCol] or next;
        # normalize bin counts for sample depth
        my @normBinCounts =
            map {roundCount($$bin[$_ + $nBinLeadCol] / $weights[$_])} 0..$maxSampleI;   
        push @$bin, @normBinCounts;         
        # the most typical sample count for this bin
        $$bin[$medCovCol] = roundCount(median(@normBinCounts)) || 1;
        # adjust bin size in case a CNV is in one sample
        $$bin[$adjLenCol] = int($$bin[$noGapLenCol] / ($$bin[$medCovCol] / $binCount) + 0.5); 
    }           
}

1;
