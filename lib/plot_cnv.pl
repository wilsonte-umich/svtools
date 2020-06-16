use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $tmpH 
            $inputDir $tmpDir $outputDir $plotDir
            %pPrpCols
            %binFileCols $nbinFileLeadCol %binColOffsets %gntCols);
my %c = %pPrpCols;
my $bsCol = $binFileCols{start};
my $beCol = $binFileCols{end};
my $cnCol = $binFileCols{cnRaw};
my ($group, $callBinsFile, $idxAnomFile, $genSnpsDir,
    $centerWidth,
    $gapFile, $segDupFile, $maxSnpCount,
    $refName, $altName,
    $widthPadding, $bamGlobXX, $minQualXX);
my ($callGroup, $callSmpFile, %samples, $nSamples, @samples, %allSamples);
my ($posH, $datH, $spnH, $snpH, $prpH, $lstH, $bamH, $gapH);

# manage options
sub setOptions_plot_cnv {
    setOptionValue(\$callBinsFile,  'call-bins-file');
    setOptionValue(\$idxAnomFile,   'idx-anom-file');
    setOptionValue(\$genSnpsDir,    'gen-snps-dir');
    setOptionValue(\$group,         'group');
    setOptionValue(\$plotDir,       'plot-dir');
    setOptionValue(\$tmpDir,        'tmp-dir',      '/tmp');
    setOptionValue(\$centerWidth,   'center-width', 5);
    setOptionValue(\$gapFile,       'gap-file');
    setOptionValue(\$segDupFile,    'seg-dup-file');
    setOptionValue(\$maxSnpCount,   'max-snp-count', 5);
    setOptionValue(\$refName,       'ref-name', "Ref");
    setOptionValue(\$altName,       'alt-name', "Alt");
    #------------------------------------------------------
    -e $callBinsFile or die "$error: $callBinsFile not found\n";
    $callBinsFile =~ m/svtools\.call\.bins\.(.+)\.bgz$/ or die "not a call bins file: $callBinsFile\n";
    $callGroup = $1;
    $callSmpFile = $callBinsFile;
    $callSmpFile =~ s/bins\.$callGroup\.bgz/samples\.$callGroup\.txt/;
    -e $idxAnomFile or die "$error: $idxAnomFile not found\n";
    $idxAnomFile =~ m/svtools\.index\.anomalies\.(.+)\.bgz$/ or die "not an index anomalies file: $idxAnomFile\n";
    -d $plotDir or die "$error: $plotDir does not exist or is not a directory\n";
    -d $tmpDir  or die "$error: $tmpDir does not exist or is not a directory\n";
    -d $genSnpsDir or die "$error: $genSnpsDir does not exist or is not a directory\n";
    $centerWidth =~ m|\D| and die "$error: center-width must be an integer number\n";
    $gapFile and (-e $gapFile or die "$error: $gapFile not found\n");
    $segDupFile and (-e $segDupFile or die "$error: $segDupFile not found\n");
    $outputDir = $plotDir;
    $widthPadding = ($centerWidth - 1) / 2;
}

# main execution block
sub svtools_plot_cnv {

    # initialize
    setOptions_plot_cnv();    
    print STDERR "$utility plot_cnv: " . getTime(), "\n";
    openInputStream($callSmpFile, \$tmpH);
    my $samples = <$tmpH>;
    closeHandles($tmpH);
    chomp $samples;
    @samples = split("\t", $samples);
    %samples = map { $samples[$_] => $_ } 0..$#samples;
    $nSamples = @samples;

    # collect data for cnv span on the input stream
    print STDERR "collecting data\n";
    my $plotFile = getTmpFile("plot_cnv", 'data', "gz");
    openOutputStream($plotFile, \$datH, $TRUE);
    my $spanFile = getTmpFile("plot_cnv", 'spans', "gz");
    openOutputStream($spanFile, \$spnH, $TRUE);
    my $snpFile = getTmpFile("plot_cnv", 'snps', "gz");
    openOutputStream($snpFile, \$snpH, $TRUE);
    my $gapDupFile  = getTmpFile("plot_cnv", 'gaps_dups', "gz");
    openOutputStream($gapDupFile, \$gapH, $TRUE);
    my $lstFile = getOutFile("list", $group, "txt");
    openOutputStream($lstFile, \$lstH, $FALSE);
    openInputStream(undef, \$prpH);    
    my $header = <$prpH>;
    print $lstH $header;
    my $nCnvs = 0;
    while (my $prp = <$prpH>) {
        $nCnvs++;
        print $lstH $prp;
        chomp $prp;
        my @prp = split("\t", $prp);
        my $width = $prp[$c{end}] - $prp[$c{start}] - 1;
        my $start = $prp[$c{start}] - $width * $widthPadding;
        my $end = $prp[$c{end}] + $width * $widthPadding;
        my $plotRegion = "$prp[$c{chrom}]:$start-$end";   
        $prp[$c{region}] =~ s/:/-/;
        print STDERR "  $prp[$c{grpId}]  = $plotRegion\n";
        printDepthBins(\@prp, $plotRegion, $start, $end, $prp[$c{start}], $prp[$c{end}]);
        printAnomSpans(\@prp, $plotRegion, $start, $end);
        printSNPs(\@prp, $plotRegion);
        printGaps(\@prp, $plotRegion);
    }
    closeHandles($datH, $spnH, $snpH, $prpH, $lstH, $gapH);

    # make the plots
    if ($nCnvs) {
        print STDERR "making plots\n";
        $ENV{plotFile}     = $plotFile;
        $ENV{spanFile}     = $spanFile;
        $ENV{snpFile}      = $snpFile;
        $ENV{gapDupFile}   = $gapDupFile;
        $ENV{plotDir}      = $plotDir;
        $ENV{centerWidth}  = $centerWidth;
        $ENV{maxSnpCount}  = $maxSnpCount;
        $ENV{refName}      = $refName;
        $ENV{altName}      = $altName;
        system("Rscript $libDir/plot_cnv.R");
    } else {
        print STDERR "nothing to plot\n";  
    }

    # finish up 
    unlink $plotFile;
    unlink $spanFile;
    unlink $snpFile;
    unlink $gapDupFile;
}

# get the data for the coverage depth traces
sub printDepthBins {
    my ($prp, $plotRegion, $start, $end, $left, $right) = @_;
    open my $tbxH, "-|", "tabix $callBinsFile $plotRegion" or
        die "could not open tabix stream: tabix $callBinsFile $plotRegion\n";
    while (my $bin = <$tbxH>) {
        chomp $bin;
        my @bin = split("\t", $bin);
        my $center = int($bin[$bsCol] + ($bin[$beCol] - $bin[$bsCol]) / 2 + 0.5);
        my $modCN  = $bin[$cnCol];
        foreach my $sampleI(0..$#samples){
            my $dCN = $bin[$nbinFileLeadCol + $sampleI + $nSamples * $binColOffsets{smpDCNRaw}];
            print $datH join("\t",
                (map {$$prp[$c{$_}]} qw(chrom region group grpId cnvType cnvSize cnvFrac sample othSample)),
                $start, $end, $left, $right,
                $samples[$sampleI], $center, $modCN, $dCN), "\n";            
        }        
    }
    close $tbxH;
}

# get the data for the anomalous fragment tracks
sub printAnomSpans {
    my ($prp, $plotRegion, $start, $end) = @_;
    
    # split by sample and pair type      
    open my $tbxH, "-|", "tabix $idxAnomFile $plotRegion" or
        die "could not open tabix stream: tabix $idxAnomFile $plotRegion\n";
    my %pairs;
    my @pairTypes = qw(del dup invF invR);    
    my %pairTypes = map { $_ => 1 } @pairTypes;
    %allSamples = map { $_ => 1 } @samples;
    while (my $bed = <$tbxH>) {
        chomp $bed;
        my @bed = split("\t", $bed);
        my $okSize = (($bed[2] - $bed[1]) < 2.5e6);
        my $inWindow = (($bed[1] >= $start and $bed[1] <= $end) and 
                        ($bed[2] >= $start and $bed[2] <= $end));
        my $oneEndIn = (($bed[1] >= $start and $bed[1] <= $end) or 
                        ($bed[2] >= $start and $bed[2] <= $end));
        ($inWindow or ($oneEndIn and $okSize)) or next;     
        $bed[3] =~ m/(.+)\..+,(\w+)/; # sample, cnvType        
        $pairTypes{$2} or next;
        push @{$pairs{$1}{$2}{$bed[1]}{$bed[2]}}, \@bed;
        $allSamples{$1}++;
    }    
    close $tbxH;
    
    # plot a track for each pair type
    foreach my $sample(keys %allSamples){
        my $maxLineAll = 1;
        foreach my $pairType(@pairTypes){
            my $bed = $pairs{$sample}{$pairType} or next;
            my ($pckBed, $maxLine) = packSpans($bed);
            foreach my $bed(@$pckBed){
                print $spnH join("\t",
                    (map {$$prp[$c{$_}]} qw(grpId sample)),
                    $sample, $pairType, @$bed[2..3], $$bed[0]+$maxLineAll
                ),"\n"; 
            }
            $maxLineAll += $maxLine + 2;
        }    
    }    
}
sub packSpans {  # minimize vertical display by packing spans
    my ($spns) = @_; # expects $$spns{$p1}{$p2} = [ [], ... ]
    my ($baseLine, $maxLine, %lineMaxs, @pckFeats) = (0, 0);
    foreach my $p1(sort {$a <=> $b} keys %$spns){
        foreach my $p2(sort {$a <=> $b} keys %{$$spns{$p1}}){
            foreach my $feat(@{$$spns{$p1}{$p2}}){
                my $line = $baseLine;
                while(defined $lineMaxs{$line} and $lineMaxs{$line} >= $p1){ $line++ }
                $maxLine >= $line or $maxLine = $line;
                $lineMaxs{$line} = $p2;
                unshift @$feat, $line; # put assigned line for this feature in array index 0
                push @pckFeats, $feat; # assemble an array of packed features
            }  
        }   
    }
    return (\@pckFeats, $maxLine); # return the packed features and maximum assigned line
}

# print SNP coverage
sub printSNPs {
    my ($prp, $plotRegion) = @_;    
    foreach my $sample(keys %allSamples){
    #foreach my $sample($$prp[$c{sample}]){   
        my $genSnpsFile = "$genSnpsDir/svtools.genotype.snps.$sample.bgz";
        -e $genSnpsFile or next;
        open my $tbxH, "-|", "tabix $genSnpsFile $plotRegion" or
            die "could not open tabix stream: tabix $genSnpsFile $plotRegion\n";      
        while (my $snp = <$tbxH>) {
            chomp $snp;
            my @snp = split("\t", $snp);
            print $snpH join("\t",
                (map {$$prp[$c{$_}]} qw(grpId sample)),
                $sample,
                (map {$snp[$gntCols{$_}]} qw (end nRef nAlt refBase altBase))
            ),"\n";
        }        
        close $tbxH;
    }         
}

# highlight genome gaps and seg dups
sub printGaps {
    my ($prp, $plotRegion) = @_;
    print $gapH join("\t", # make sure the file is not empty
        (map {$$prp[$c{$_}]} qw(grpId)),
        -1, -1, 'black'
    ), "\n"; 
    $gapFile and printGaps_($gapFile, 'darkgrey', $prp, $plotRegion);
    $segDupFile and printGaps_($segDupFile, 'lightblue', $prp, $plotRegion);
}
sub printGaps_ {
    my ($file, $col, $prp, $plotRegion) = @_;
    open my $tbxH, "-|", "tabix $file $plotRegion" or
        die "could not open tabix stream: tabix $file $plotRegion\n";    
    while (my $bed = <$tbxH>) {
        chomp $bed;
        my @bed = split("\t", $bed);
        print $gapH join("\t",
            (map {$$prp[$c{$_}]} qw(grpId)),
            $bed[1], $bed[2], $col
        ), "\n"; 
    }    
    close $tbxH;
}

1;

#Col	Field	Description
#1	QNAME	Query template/pair NAME
#2	FLAG	bitwise FLAG
#3	RNAME	Reference sequence NAME
#4	POS	1-based leftmost POSition/coordinate of clipped sequence
#5	MAPQ	MAPping Quality (Phred-scaled)
#6	CIAGR	extended CIGAR string
#7	MRNM	Mate Reference sequence NaMe (?=? if same as RNAME)
#8	MPOS	1-based Mate POSistion
#9	TLEN	inferred Template LENgth (insert size)
#10	SEQ	query SEQuence on the same strand as the reference
#11	QUAL	query QUALity (ASCII-33 gives the Phred base quality)
#12+	OPT	variable OPTional fields in the format TAG:VTYPE:VALUE
#Each bit in the FLAG field is defined as:
#
#Flag	Chr	Description
#0x001	p	the read is paired in sequencing
#0x002	P	the read is mapped in a proper pair
#0x004	u	the query sequence itself is unmapped
#0x008	U	the mate is unmapped
#0x010	r	strand of the query (1 for reverse)
#0x020	R	strand of the mate
#0x040	1	the read is the first read in a pair
#0x080	2	the read is the second read in a pair
#0x100	s	the alignment is not primary
#0x200	f	the read fails platform/vendor quality checks
#0x400	d	the read is either a PCR or an optical duplicate
#0x800      supplementary alignment
