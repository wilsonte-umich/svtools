use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir $slurp
            $TRUE $FALSE
            $tmpH $inH $outH
            $inputDir $tmpDir $outputDir $plotDir
            %gntCols);
my (@globs, $group, $maxMem);
my %c = %gntCols;
my (%counts, %fracs);

# manage options
sub setOptions_pool {
    setOptionValue(\$group,      'group');
    setOptionValue(\$inputDir,   'input-dir');
    setOptionValue(\$tmpDir,     'tmp-dir',      '/tmp');
    setOptionValue(\$maxMem,     'max-mem',      1000000000);
    setOptionValue(\$outputDir,  'output-dir');     
    #-------------------------------------
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    -d $outputDir or die "$error: directory not found: $outputDir\n";   
}

# main execution block
sub svtools_pool {

    # initialize
    setOptions_pool();    
    print STDERR "$utility pool: " . getTime(), "\n";
    
    # pool SNP counts
    my $inStream = "gunzip -c $inputDir/svtools.genotype.snps.*.bgz";
    openInputStream($inStream, \$inH);
    my %snps;
    while (my $line = <$inH>) {
        chomp $line;
        my @f = split("\t", $line);
        my $snp = join(":", @f[$c{chrom}, $c{end}, $c{snpId}, $c{refBase}, $c{altBase}]);
        $snps{$snp} or $snps{$snp} = [0, 0];
        $snps{$snp}[0] += $f[$c{nRef}];
        $snps{$snp}[1] += $f[$c{nAlt}];   
    }
	closeHandles($inH);    
    
    # report the results    
	my $poolFile = getOutFile("snps", $group, "bgz");
    openOutputStream($poolFile, \$outH, undef, undef,
            "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c");
    foreach my $snp(keys %snps){
        my ($chrom, $pos, $snpId, $refBase, $altBase) = split(":", $snp);
        my $nReads = $snps{$snp}[0] + $snps{$snp}[1];
        print $outH join("\t",
            $chrom, $pos-1, $pos, $snpId,
            $nReads, '.',
            $refBase, $altBase,
            $snps{$snp}[0], $snps{$snp}[1]
        ), "\n";
        $counts{$nReads}++;
        if ($nReads) {
            my $frac = int(($snps{$snp}[1] / $nReads) / 0.05 + 0.5) * 0.05;
            $fracs{$frac}++
        }
    }
    closeHandles($outH);
    
    # report tallies
    printGenotypeSummary(\%counts, 'coverage',    qw(coverage nSNPs));
    printGenotypeSummary(\%fracs,  'frequencies', qw(fracAlt  nSNPs));
    
    # finish up
    system("tabix -p bed $poolFile");
}

sub printGenotypeSummary {
    my ($hash, $name, @lbls) = @_;
	my $outFile = getOutFile($name, $group, "txt");
    openOutputStream($outFile, \$outH);    
    print $outH join("\t", @lbls), "\n";
    foreach my $val(sort {$a <=> $b} keys %$hash){
        print $outH join("\t", $val, $$hash{$val}), "\n";
    }
    closeHandles($outH);
}

1;
