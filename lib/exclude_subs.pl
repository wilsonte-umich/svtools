use strict;
use warnings;

#========================================================================
# 'exclude.pl' supplies the exlude function, adapted from bed_util crossing
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($utility $command);
my $error;
my (%tmpExclRegions, %exclRegions, $exclBinSize);
my ($collRefBases, $maxCollRefLength, $maxCollRefEnd, $nCollRefFeatures) = (0, 0, 0, 0);
#========================================================================


#========================================================================
# collect and process the excluded, i.e. reference, features
#------------------------------------------------------------------------
sub initializeExclude {
    my ($excludefile) = @_;
    $excludefile or return;
    $error = "$utility $command error";
    loadBedFile($excludefile, \%tmpExclRegions);
    collapseRefFeatures();
    %tmpExclRegions = ();    
    my $aveCollRefLength = int($collRefBases / $nCollRefFeatures);
    $exclBinSize = int($maxCollRefEnd / 100);  # nominally split largest chromosome into 100 bins for faster searching (but higher memory use)
    my $minRefBinSize = $aveCollRefLength * 2;
    $exclBinSize = $exclBinSize > $minRefBinSize ?  $exclBinSize : $minRefBinSize;   # but use larger bins if reference features are especially long
    binRefFeatures();       
}
#========================================================================


#========================================================================
# check whether a feature should be excluded
#------------------------------------------------------------------------
sub checkExclude {  
    my ($chrom, $testStart, $testEnd) = @_;
    if(!$exclRegions{$chrom}){
        $chrom = "chr$chrom";
        $exclRegions{$chrom} or return;
    }
    #($testEnd - $testStart) > 0 or die "$error: malformed request to checkExclude: $chrom, $testStart, $testEnd\n";
    if ($testEnd - 1 < $testStart) {
        ($testStart, $testEnd) = ($testEnd - 1, $testStart + 1);
    }
    my ($startBin, $endBin) = getCrossedBins($testStart, $testEnd); 
    for (my $bin = $startBin; $bin <= $endBin; $bin += $exclBinSize){  # check all bins crossed by query feature
        $exclRegions{$chrom}{$bin} or next;
        foreach my $refRegion(@{$exclRegions{$chrom}{$bin}}){  
            $$refRegion[0] <= $testEnd and $$refRegion[1] >= $testStart and return 1; # any overlap with exclude region excludes the test feature            
        }
    }
    return undef;
}
sub countExcluded {  
    my ($chrom, $testStart, $testEnd) = @_;
    if(!$exclRegions{$chrom}){
        $chrom = "chr$chrom";
        $exclRegions{$chrom} or return 0;
    }
    #($testEnd - $testStart) > 0 or die "$error: malformed request to checkExclude: $chrom, $testStart, $testEnd\n";
    if ($testEnd - 1 < $testStart) {
        ($testStart, $testEnd) = ($testEnd - 1, $testStart + 1);
    }
    my ($startBin, $endBin) = getCrossedBins($testStart, $testEnd);
    my $sum = 0;
    for (my $bin = $startBin; $bin <= $endBin; $bin += $exclBinSize){  # check all bins crossed by query feature
        $exclRegions{$chrom}{$bin} or next;
        foreach my $refRegion(@{$exclRegions{$chrom}{$bin}}){
            ($$refRegion[0] <= $testEnd and $$refRegion[1] >= $testStart) or next;
            if($$refRegion[0] <= $testStart and $$refRegion[1] >= $testEnd){
                $sum += $testEnd - $testStart;
            } elsif($$refRegion[0] >= $testStart and $$refRegion[1] <= $testEnd) {
                $sum += $$refRegion[1] - $$refRegion[0];
            } elsif($$refRegion[0] < $testStart) {
                $sum += $$refRegion[1] - $testStart;
            } elsif($$refRegion[1] > $testEnd) {
                $sum += $testEnd - $$refRegion[0];
            }          
        }
    }
    return $sum;
}
#========================================================================


#========================================================================
# worker subs
#------------------------------------------------------------------------
sub loadBedFile {
    my ($bedFile, $featuresHash) = @_;
    open my $inH, "<", $bedFile or die "$error: could not open $bedFile: $!\n";
    while (my $line = <$inH>){
        $line =~ m|^\s*#| and next;  # ignore comment lines
        chomp $line;
        $line =~ s/\r//g;
        my ($chrom, $start, $end) = split("\t", $line, 4);
        $chrom or next;  # ignore blank lines                
        parseInt(\$start);  # start and end must be present and integer numbers
        parseInt(\$end);               
        $$featuresHash{$chrom}{$start}{$end}++;
    }  
    close $inH;
}
sub parseInt {  # validate and uncommify integer values
    my ($int) = @_;  # int passed as reference
    defined $$int or die "$error: missing start or end in exclude-file\n";
    $$int =~ s/,//g;
    $$int =~ m|\D| and die "$error: invalid integer in exclude-file\n";
}
#------------------------------------------------------------------------
sub collapseRefFeatures {  # convert reference features to region hash, collapsing for certain score types
    foreach my $chrom(keys %tmpExclRegions){    
        my ($regionEnd, $regionStart)= (0); 
        foreach my $start(sort {$a <=> $b} keys %{$tmpExclRegions{$chrom}}){     
            if(!$regionEnd){
                $regionStart = $start; 
            } elsif($start > $regionEnd){
                push @{$exclRegions{$chrom}}, [$regionStart, $regionEnd];
                ($regionStart, $regionEnd) = ($start, 0);
            } 
            foreach my $end(keys %{$tmpExclRegions{$chrom}{$start}}){
                $regionEnd >= $end or $regionEnd = $end;
            }
        }
        push @{$exclRegions{$chrom}}, [$regionStart, $regionEnd];
        $exclRegions{$chrom} or next;   
        foreach my $refRegion(@{$exclRegions{$chrom}}){  # collect reference feature stats
            my ($start, $end) = @$refRegion;
            my $featureLength = $end - $start;
            $collRefBases += $featureLength;
            $maxCollRefLength >= $featureLength or $maxCollRefLength = $featureLength;
            $maxCollRefEnd >= $end or $maxCollRefEnd = $end;
            $nCollRefFeatures++;
        }         
    }     
}
#------------------------------------------------------------------------
sub binRefFeatures {  # split reference features into bins on each chromosome strand to speed up crossing search
    my %binned;
    foreach my $chrom(keys %exclRegions){   
        foreach my $refRegion(@{$exclRegions{$chrom}}){
            my ($startBin, $endBin) = getCrossedBins(@$refRegion);
            for (my $bin = $startBin; $bin <= $endBin; $bin += $exclBinSize){
                push @{$binned{$chrom}{$bin}}, $refRegion;
            }
        }
    } 
    %exclRegions = %binned;
}
sub getCrossedBins {  # the lowest and highest bins crossed by a region
    my ($start, $end) = @_;
    return (int($start     / $exclBinSize) * $exclBinSize,
            int(($end - 1) / $exclBinSize) * $exclBinSize);  # end bins are converted to 0-referenced start coordinates
}
#========================================================================

1;
