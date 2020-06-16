use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $rptH
            $inputDir $tmpDir $outputDir $plotDir
            $_gap $_split $_clip $pTLenBin);
my (@globs, $group, $delete, $maxMem, $chroms);
my ($gZip, $stepsAfter, $stepsBefore);
my (%chroms, @chroms, @samples);

# manage options
sub setOptions_cat {
    setOptionValue(\$group,     'group');
    setOptionValue(\$inputDir,  'input-dir');
    setOptionValue(\$delete,    'delete');
    setOptionValue(\$outputDir, 'output-dir');
    setOptionValue(\$tmpDir,    'tmp-dir',         '/tmp');
    setOptionValue(\$maxMem,    'max-mem',         1000000000);
    setOptionValue(\$chroms,    'chroms');
    #-------------------------------------
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    $maxMem  =~ m|\D| and die "$error: option --max-mem must be an integer number of RAM bytes\n";    
    %chroms = map { $_ => 1} split(",", $chroms);
    @chroms = sort {$a cmp $b} keys %chroms;
}

# main execution block
sub svtools_cat {
    # take no files as input, collected based on group and inputDir

    # initialize streams
    setOptions_cat();    
    print STDERR "$utility cat: " . getTime(), "\n";
    my $smpFile = getOutFile("samples",   $group, "txt");
    my $setFile = getOutFile("sets",      $group, "bgz");
    my $anmFile = getOutFile("anomalies", $group, "bgz");
    openOutputStream($smpFile, \my $smpH);
    my $bgZip = "sort -S $maxMem"."b -k1,1 -k2,2n -k3,3n | bgzip -c";
    openOutputStream($setFile, \my $setH, undef, undef, $bgZip);
    openOutputStream($anmFile, \my $anmH, undef, undef, $bgZip);
    my $nCNVLines = catFiles('sets', $setH);
    my $nAnmLines = catFiles('anomalies', $anmH);
    print $smpH join("\t", @samples), "\n";
    closeHandles($smpH, $setH, $anmH);
    system("tabix -p bed $setFile");
    system("tabix -p bed $anmFile");
    my $nCNVs = $nCNVLines / 2;
    my $nAnms = $nAnmLines;
    print STDERR "$utility cat: $nAnms anomalies in $nCNVs CNVs for group $group\n";
}

sub catFiles {
    my ($type, $outH) = @_;
    my $nLines = 0;
    foreach my $inFile(glob(getInFile("find", $type, "$group.*", "gz"))){
        $inFile =~ m/svtools\.find\.$type\.$group\.(\w+)/ or
            die "$utility cat error: $inFile is not a valid $type file name for $group\n";    
        my $chrom = $1;
        $chroms{$chrom} or next;
        openInputStream($inFile, \$inH, $FALSE, $FALSE, $TRUE);
        my $line = <$inH>;
        chomp $line;
        if($type eq 'sets'){
            $line =~ m/^#(.+)/ or die "$utility cat error: no header line in bin file: $inFile\n";
            my %tags = getFileHeader($1, 'find', qw(group samples));
            @samples = split(",", $tags{samples});            
        }
        while (my $line = <$inH>){
            print $outH $line;
            $nLines++;
        } 
        closeHandles($inH);
    }
    return $nLines;
}

1;
