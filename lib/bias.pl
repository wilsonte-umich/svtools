use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir $slurp
            $TRUE $FALSE
            $tmpH $inH
            $inputDir $tmpDir $outputDir $plotDir
            %gntCols);
my %c = %gntCols;
my $stateCol = $c{readBases};
my (@globs, $sample, $snpFile, $maxMem, $chroms,
	$zeroProb, $persistence, $errorRate);
my ($biasH, $cnvH) = @_;
my (%chroms, @snps, %snpCov, $prevF);
my ($smpSumRef, $smpSumAlt, $smpSumN, $cnvId) = (0, 0, 0, 1);
my @biasTypes = qw(het refDup altDup refDel altDel);

# manage options
sub setOptions_bias {
    setOptionValue(\$sample,     'sample');
    setOptionValue(\$inputDir,   'input-dir');
    setOptionValue(\$snpFile,    'snp-file');
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);
    setOptionValue(\$outputDir,  'output-dir');
    setOptionValue(\$chroms,     'chroms');
	setOptionValue(\$zeroProb,   'zero-prob',		0.99);
	setOptionValue(\$persistence,'persistence',		0.995);
	setOptionValue(\$errorRate,  'base-error-rate', 0.001);	
    #-------------------------------------
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    -e $snpFile or die "$error: $snpFile does not exist\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";
    %chroms = map { $_ => 1} split(",", $chroms);
}

# main execution block
sub svtools_bias {

    # initialize
    setOptions_bias();    
    print STDERR "$utility bias: " . getTime(), "\n";

    # open input and output SNP streams
    print STDERR "finding runs of SNP genotype bias\n";
    my $gntFile = getInFile("genotype", "snps", $sample, "bgz");
    openInputStream("gunzip -c $gntFile", \$inH);
	my $biasFile = getOutFile("snps", $sample, "bgz");
    openOutputStream($biasFile, \$biasH, undef, undef,
            "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
	my $cnvFile = getOutFile("cnvs", $sample, "bgz");
    openOutputStream($cnvFile, \$cnvH, undef, undef,
            "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    
    # run SNP bias analysis on each chromosome in turn    
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        $chroms{$f[$c{chrom}]} or next;
        $prevF and $$prevF[$c{chrom}] ne $f[$c{chrom}] and runBiasHMM();
        $smpSumN   += $f[$c{nReads}];
        $smpSumAlt += $f[$c{nAlt}];
        push @snps, \@f;
        $prevF = \@f;
    }
    runBiasHMM();
    closeHandles($inH, $biasH, $cnvH);

    # finish up
    system("tabix -p bed $biasFile");
    system("tabix -p bed $cnvFile");
}

sub runBiasHMM {
    @snps or next;
    print STDERR "  $$prevF[$c{chrom}]\n";
    
    # report the whole-chromosome fracAlt, which need not == 0.5 if there is a mapping bias
    my $smpFracAlt = $smpSumAlt / $smpSumN;
	print STDERR "    $smpFracAlt = fraction of alternative alleles\n";
    ($smpFracAlt < 0.1 or $smpFracAlt > 0.9) and return printHemizygous();
    $smpSumRef = $smpSumN - $smpSumAlt;

	# report lambda, i.e. the average number of reads per snp
	my $nSnps = @snps;
	my $lambda = $smpSumN / $nSnps;
	print STDERR "    $nSnps = # of SNPs\n";	
	print STDERR "    $smpSumN = # of crossing reads\n";
	print STDERR "    $lambda = mean reads per SNP\n";	

    # print the chromosomes snps to a temporary file for R
    my ($refSumRef, $refSumAlt, $refSumN) = loadSnpCounts();
    my $tmpFile = getTmpFile('bias_HMM', $sample, "gz");
    openOutputStream($tmpFile, \my $tmpH, $TRUE);
    foreach my $i(0..$#snps){
        my $refCov  = $snpCov{$snps[$i][$c{end}]};
        print $tmpH join("\t",
            0, $i,
            $snps[$i][$c{nRef}], $snps[$i][$c{nAlt}],
            ($$refCov[0] + $$refCov[1]) * $smpSumN / $refSumN, # SNP-specific expected coverage (lambda)
            $$refCov[1] / ($$refCov[0] + $$refCov[1])), "\n";  # SNP-specific expected fracAlt
    }
    closeHandles($tmpH);
    
    # set the emissions probablities and run HMM segmentation
    $ENV{dataFile} = $tmpFile;
	#$ENV{lambda}   = $lambda;	
    #$ENV{fracAlt}  = $fracAlt;
	$ENV{errorRate}  = $errorRate;
    open my $segH, "-|", "Rscript $libDir/bias.R | segment -z $zeroProb -p $persistence -w"
        or die "could not open segment stream: $!\n";
	my ($sumLogProb0, $sumLogProbCnv, $runState, $runStartI) = (0, 0);
    while (my $line = <$segH>) {
        chomp $line;
        my ($state, $logProbs, $i) = split("\t", $line);
		my @logProbs = split(",", $logProbs);		
        $snps[$i][$stateCol]   = $biasTypes[$state];
		$snps[$i][$stateCol+1] = -$logProbs[0] - -$logProbs[$state];	
        print $biasH join("\t", @{$snps[$i]}), "\n";
        if(defined $runState and $runState != $state){
			$runState and printBiasCnv($runStartI, $i-1, $runState, $sumLogProb0, $sumLogProbCnv);
			$runStartI = undef;
			($sumLogProb0, $sumLogProbCnv) = (0, 0);
		}		
		$sumLogProb0   += $logProbs[0];
		$sumLogProbCnv += $logProbs[$state];	
        $runState = $state;
        !(defined $runStartI) and $runStartI = $i;
    }
	$runState and printBiasCnv($runStartI, $nSnps-1, $runState, $sumLogProb0, $sumLogProbCnv);
    closeHandles($segH);

    # finish the chromosome
    unlink $tmpFile;
    @snps = ();
    ($smpSumN, $smpSumAlt) = (0, 0);
}
sub printBiasCnv {
    my ($startI, $endI, $state, $sumLogProb0, $sumLogProbCnv) = @_;
    my ($nRef, $nAlt) = (0, 0);
    foreach my $i($startI..$endI){
        $nRef += $snps[$i][$c{nRef}];
        $nAlt += $snps[$i][$c{nAlt}];
    }
    print $cnvH join("\t",
        $$prevF[$c{chrom}], $snps[$startI][$c{start}], $snps[$endI][$c{end}],
        $cnvId, $endI-$startI+1, '.',
        $sample, $nRef, $nAlt, $biasTypes[$state], -$sumLogProb0 - -$sumLogProbCnv
    ), "\n";
    $cnvId++;     
}

# print a hemizygous/UPD chromosome 
sub printHemizygous {
    print STDERR "skipping HMM, $sample $$prevF[$c{chrom}] shows hemizygosity or uniparental disomy";
    foreach my $i(0..$#snps){
        $snps[$i][$stateCol]   = 'loh';
		$snps[$i][$stateCol+1] = 999;        
        print $biasH join("\t", @{$snps[$i]}), "\n";        
    }
}

# load reference SNP counts for current chromosome
sub loadSnpCounts {
    %snpCov = ();
    openInputStream("tabix $snpFile $$prevF[$c{chrom}]", \my $inH);
    my ($sumRef, $sumAlt) = (0, 0);
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        $snpCov{$f[2]} = [$f[8], $f[9]];
        $sumRef += $f[8];
        $sumAlt += $f[9];
    }
    closeHandles($inH);
    return ($sumRef, $sumAlt, $sumRef + $sumAlt);
}

1;
