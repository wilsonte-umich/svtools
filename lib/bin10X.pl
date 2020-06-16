use strict;
use warnings;

# define variables
use vars qw($version $utility $error $libDir
            $TRUE $FALSE
            $inH $tmpH $outH
            $inputDir $tmpDir $outputDir $plotDir);
#require "$libDir/exclude_subs.pl";
my (@globs, $sample, $chrom,
    $weightPerCell, $gapFile, $excludeFile,
    $minQual, $slurpSize);
my ($wrkCount, $maxEnd,
    $prevBinEnd, $brkStart, $brkEnd) = (0) x 20;
my $ALL = 'all_cells';
my ($nCells, @allCells, %cells) = (0, $ALL);
my (@frags, $binWeight, %nextBinCount) = ();
my ($nFragments, $nBins, @binSizes, @binCounts) = (0) x 2;
my $f = 0x001 + # the read is paired in sequencing
        0x002;  # the read is mapped in a proper pair
my $F = 0x004 + # the query sequence itself is unmapped
        0x008 + # the mate is unmapped
        0x010 + # strand of the query (1 for reverse)
        0x100 + # the alignment is not primary
        0x200 + # the read fails platform/vendor quality checks
        0x400 + # the read is either a PCR or an optical duplicate
        0x800;  # supplementary alignment
            
# manage options
sub setOptions_bin10X {
    setOptionValue(\$inputDir, 'input-dir');
    setOptionValue(\$chrom,     'chrom');
    setOptionValue(\$weightPerCell, 'weight-per-cell', 10);     
    setOptionValue(\$minQual,   'min-qual',     30);
    setOptionValue(\$slurpSize,   'buffer-size',     "2500M");
    #setOptionValue(\$gapFile,    'gap-file');
    #setOptionValue(\$excludeFile,'exclude-file');    
    setOptionValue(\$outputDir, 'output-dir');    
    setOptionValue(\$sample,    'sample');    
    setOptionValue(\$tmpDir,    'tmp-dir',     '/tmp');     
    #-------------------------------------
    -d $inputDir or die "$error: $inputDir does not exist or is not a directory\n";
    $weightPerCell =~ m|\D| and die "$error: weight-per-cell must be an integer\n";
    $minQual =~ m|\D| and die "$error: min-qual must be an integer\n";
    $excludeFile and (-e $excludeFile or die "$error: file not found: $excludeFile\n");
    -d $outputDir or die "$error: $outputDir does not exist or is not a directory\n";
    -d $tmpDir or die "$error: $tmpDir does not exist or is not a directory\n";
}

# main execution block
sub svtools_bin10X {
    (@globs) = @_;
    my $usage = "$error: bin10X does not accept input on command line; --input-dir should point to a CellRanger outs directory";
    $globs[0] and die "$error: $usage\n";

    # initialize
    setOptions_bin10X();    
    print STDERR "$utility bin10X: " . getTime(), "\n";
    #initializeExclude($excludeFile);
    
    ###############################
    #my $checkFile = getOutFile("bins", "$sample.$chrom", "gz");
    #-e $checkFile and exit;

    # load cell information from discover10X
    loadCells10X();
    
    # parse and count the bins
    # flags recover forward read of a proper pair, not secondary, supplemental or duplicate    
    my $outFile = getOutFile("bins", "$sample.$chrom", "gz");
    openOutputStream($outFile, \$outH, $TRUE);
    print $outH join("\t", qw(chrom start end), @allCells), "\n";
    run_bin10X();
    closeHandles($inH, $outH);

    # report counts
    reportCount($nFragments, "fragments processed over $nCells allowed cell barcodes"); 
    reportCount($nBins,      "bins reported");
    reportCount(roundCount(median(@binSizes), 1),      "median bin size");
    reportCount(roundCount(median(@binCounts), 1),      "median bin count over all samples");
    reportCount($binWeight,      "target bin count over all samples");
}

## load cells previously designated as "true" by discover10X
#sub loadCells10X {
#    my $cellFile = getInFile('discover10X', "cells", "$sample", "txt");
#    openInputStream($cellFile, \$inH);
#    my $header = <$inH>;
#    while(my $line = <$inH>){
#        my ($CB, $trueCell) = split("\t", $line);
#        $trueCell or next;
#        $cells{$CB}++;
#        $nCells++;
#        push @allCells, $CB;
#    }
#    closeHandles($inH);
#    $binWeight = $nCells * $weightPerCell;
#}

# load cells previously designated as "true" by CellRanger
sub loadCells10X {
    my $cellsFile = "$inputDir/per_cell_summary_metrics.csv";
    open my $inH, "<", $cellsFile or die "$error: could not open $cellsFile for reading: $!\n";
    my $header = <$inH>;
    while(my $line = <$inH>){
        my ($CB, $cellId) = split(",", $line);
        $cells{$CB} = $cellId;
        $nCells++;
        #push @allCells, $CB;
        push @allCells, $cellId;
    }
    close $inH;
    $binWeight = $nCells * $weightPerCell;
}
#barcode,cell_id,total_num_reads,num_unmapped_reads,num_lowmapq_reads,num_duplicate_reads,num_mapped_dedup_reads,frac_mapped_duplicates,effective_depth_of_coverage,effective_reads_per_1Mbp,raw_mapd,normalized_mapd,raw_dimapd,normalized_dimapd,mean_ploidy,ploidy_confidence,is_high_dimapd,is_noisy,est_cnv_resolution_mb
#AAACCTGCACCACACG-1,0,1401732,6661,222057,163321,1009693,0.11651371303501667,0.0700481200042838,370,0.11524713366463572,0.11524713366463572,1.0562023782600596,1.0562023782600596,1.949662924883774,-2,0,0,0.9989746093750012

# parse and count the bins
sub run_bin10X {
    print STDERR "$utility bin10X: parsing and counting bins on chrom $chrom\n";
    my $bamFile = "$inputDir/possorted_bam.bam";
    openBamStream($bamFile, \$inH, $slurpSize, $minQual, $f, $F, $chrom); 
    while (my $line = <$inH>){
        # finish filtering to countable reads
        my @f = split("\t", $line, 10);
        $f[8] > 0 or next; # this is the leftmost read of a pair (i.e F/R)
        
        # only allow true cells to continue
        $f[9] =~ m/CB:Z:(\S+)/ or next; # the 10X error corrected, validate cell barcode
        defined $cells{$1} or next;
        $f[0] = $cells{$1};        
        $nFragments++;              

        # parse additional fragment information
        $f[1] = $f[3] + $f[8] - 1; # the end coordinate of the fragment
        $maxEnd >= $f[1] or $maxEnd = $f[1];    
        $wrkCount++;
        
        # as needed, handle variable-width/equal-weight bin processing
        if($wrkCount >= $binWeight){ # this fragment is guaranteed to create a break condition
            
            # this is the _first_ fragment that will create a break condition (index fragment)
            if(!$brkEnd){ 
                $brkStart = $f[3]; # break cannot occur before this fragment starts           
                $brkEnd = $maxEnd; # but could break as far as the end of any fragment until now
                
            # the first fragment that CANNOT contribute to the previous bin
            } elsif($f[3] > $brkEnd){
                
                # initialize counts
                my %binCount = %nextBinCount; # the count carried over from the last bin's break
                %nextBinCount = map { $_ => 0 } @allCells; # the count to carry over to the next bin
                my ($brkCount, $binEnd, %brkSpans) = (0, 0);
                
                # examine every fragment that might contribute to bin
                foreach my $x(@frags){
                    
                    # these fragments are entirely within this bin, before the break
                    if($$x[1] < $brkStart){ 
                        $binCount{$ALL}++; # just count them
                        $binCount{$$x[0]}++;
                        
                    # these fragments are likely to be split by the break    
                    } else { 
                        my $inc = 1 / $$x[8];
                        my $CB = $$x[0];
                        my ($minPos, $maxPos, $brkWeight);
                        
                        # some fragments must have been encountered BEFORE the index fragment 
                        if($$x[3] < $brkStart){ 
                            my $bin = ($brkStart - $$x[3]) * $inc; 
                            $binCount{$ALL} += $bin;
                            $binCount{$CB} += $bin; # these bases always in this bin
                            $brkCount += 1 - $bin;
                            ($minPos, $maxPos, $brkWeight) = ($brkStart, $$x[1], 1 - $bin); # these bases in break region
                        
                        # some fragments must have been encountered AFTER the index fragment
                        } elsif($$x[1] > $brkEnd){ 
                            my $next = ($$x[1] - $brkEnd) * $inc;
                            $nextBinCount{$ALL} += $next; # these bases always in next bin
                            $nextBinCount{$CB} += $next;
                            $brkCount += 1 - $next;
                            ($minPos, $maxPos, $brkWeight) = ($$x[3], $brkEnd, 1 - $next); # these bases in break region

                        # some fragments may reside entirely in the break region
                        } else {
                            $brkCount++;
                            ($minPos, $maxPos, $brkWeight) = ($$x[3], $$x[1], 1);
                        }
                        push @{$brkSpans{$CB}}, [$minPos, $maxPos, $brkWeight];                        
                    }
                }
                
                # interpolate the break position based on actual and needed break region counts
                # this provides an estimate of the true break position, but its good enough
                # the error results in a slight variance in total bin counts
                # but each sample count in the final bin is accurate
                $binEnd = roundCount($brkStart + ($brkEnd - $brkStart) * ($binWeight - $binCount{$ALL}) / $brkCount, 1);
                $binEnd - $prevBinEnd < 0 and die "$error: parsing resulted in negative bin size\n";
                
                # add the portions of fragments within the break region to this bin or next bin
                foreach my $CB(keys %brkSpans){
                    foreach my $brkSpan(@{$brkSpans{$CB}}){
                        my ($minPos, $maxPos, $brkWeight) = @$brkSpan;
                        if ($maxPos <= $binEnd){
                            $binCount{$ALL} += $brkWeight;
                            $binCount{$CB} += $brkWeight;
                        } elsif($minPos > $binEnd){
                            $nextBinCount{$ALL} += $brkWeight;
                            $nextBinCount{$CB} += $brkWeight;
                        } else {
                            my $bin = $brkWeight * ($binEnd - $minPos) / ($maxPos - $minPos);
                            $binCount{$ALL} += $bin;
                            $binCount{$CB} += $bin;
                            $nextBinCount{$ALL} += 1- $bin;
                            $nextBinCount{$CB} += 1- $bin;
                        }   
                    }
                }

                # commit the newly finished bin
                print $outH join("\t", $chrom, $prevBinEnd, $binEnd);
                foreach my $CB(@allCells){
                    print $outH "\t".roundCount2($binCount{$CB} || 0);                    
                }
                print $outH "\n";
                push @binSizes, $binEnd - $prevBinEnd;
                push @binCounts, $binCount{$ALL} || 0;
                $prevBinEnd = $binEnd;                
                
                # handle rare regions with extremely high read density
                # "bad" genome regions that have not been excluded by user
                while($nextBinCount{$ALL} >= $binWeight){
                    my $badStart = min($binEnd + 1, $brkEnd);
                    my $badFrac = $binWeight / $nextBinCount{$ALL};
                    $binEnd = roundCount($badStart + ($brkEnd - $badStart) * $badFrac, 1);
                    $binEnd - $prevBinEnd < 0 and die "$error: parsing resulted in negative bin size while fixing a bad region\n";
                    print $outH join("\t", $chrom, $prevBinEnd, $binEnd);
                    foreach my $CB(@allCells){
                        my $count = ($nextBinCount{$CB} || 0) * $badFrac;
                        print $outH "\t".roundCount2($count);
                        $nextBinCount{$CB} -= $count; 
                    }
                    print $outH "\n";
                    $prevBinEnd = $binEnd;
                }

                # reset for the next bin
                ($wrkCount, $maxEnd, $brkEnd) = ($nextBinCount{$ALL}, $f[1], 0);
                @frags = ();
                
                # provide progress feedback
                $nBins++;
                if($nBins % 100 == 0){
                    my $b = commify($nBins);
                    my $f = commify($nFragments);
                    my $e = commify($prevBinEnd);
                    print STDERR "\t$b bins / $f fragments / $e bp\n";
                } 
            }
        }
        
        # collect the set of working fragments
        push @frags, \@f;    
    } 
}

1;

#Col	Field	Description
#1	QNAME	Query template/pair NAME
#2	FLAG	bitwise FLAG
#3	RNAME	Reference sequence NAME
#4	POS	1-based leftmost POSition/coordinate of clipped sequence
#5	MAPQ	MAPping Quality (Phred-scaled)
#6	CIAGR	extended CIGAR string
#7	MRNM	Mate Reference sequence NaMe (= if same as RNAME)
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

# both reads of a pair are marked as duplicates
#A00228:250:H3TVLDSXX:1:1110:25138:24925 1024,1024
#A00228:250:H3TVLDSXX:1:1139:28248:11710 1024,1024
#A00228:250:H3TVLDSXX:1:1146:8477:8469   1024,1024
#A00228:250:H3TVLDSXX:1:1151:9643:32659  1024,1024

# the two reads for one of those pairs
#A00228:250:H3TVLDSXX:1:1139:28248:11710 1123    22      16057320        27      84M     =       16057594        374    GTAAAATGTAGATTTTGGCACACGGGGACCTAAGACAATCCTGCCTGAAACAGGGCCTGGCATCCCATCAGTGCTCAATAAACA     FFFFF:F:FFFFFFFFF:FFFF,,FF:FFFFFFFFFF:FF:FFF,FFFFF:FFFFF,FFF,FFF,FF,,FF,FFFFFFFFFFFF    NM:i:0  MD:Z:84 AS:i:84 XS:i:84 XA:Z:14,-19785603,84M,0;       CR:Z:AAGCCGCGTTGACGGA    CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:ACCCTCCT   QT:Z:,F::FFF:   GP:i:698302392 MP:i:698302766   RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#A00228:250:H3TVLDSXX:1:1139:28248:11710 1171    22      16057594        27      100M    =       16057320        -374   TGAGAAGCTGTCCTGGTCTAAGTATCTTTCTCCCATTTTACATAAAGGAATACACAGTGTCAGAAGGAGGACCTGTGTCCAGCCCCTGTGTTCCCCACTT     FFFFFFFF,FF:FFFF:FFFFFF:FFFFF:FFFFFF:,,,FF,::FF:FF,FFFFFF:FF:,FF:FFFFFFF:FFF,:FFFFF:FF,FFF:FFFFFFF:F    NM:i:1  MD:Z:99C0       AS:i:99XS:i:94  XA:Z:14,+19785313,100M,2;       CR:Z:AAGCCGCGTTGACGGA   CY:Z:FFFFFFFFFFFFFFFF   CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:ACCCTCCT   QT:Z:,F::FFF:   GP:i:698302766  MP:i:698302392  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1

# three different reads with the same barcode and start position
# 1123 = 99 | 1024
# thus ONE of the duplicate is NOT marked as duplicate (as should be, but documentation was badly worded)
#A00228:250:H3TVLDSXX:1:1139:27434:9925  99      22      16057320        27      84M     =       16057594        374    GTAAAATGTAGATTTTGGCACACGGGGACCTAAGACAATCCTGCCTGAAACAGGGCCTGGCATCCCATCAGTGCTCAATAAACA     FFF:FFFFF:FFF:FFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFF,FFFF:F,,F,:FF,
#NM:i:0  MD:Z:84 AS:i:84 XS:i:84 XA:Z:14,-19785603,84M,0;       CR:Z:AAGCCGCGTTGACGGA    CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:ACCCTCCT   QT:Z:FFFF,F,,   GP:i:698302392 MP:i:698302766   RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1

#A00228:250:H3TVLDSXX:1:1139:28248:11710 1123    22      16057320        27      84M     =       16057594        374    GTAAAATGTAGATTTTGGCACACGGGGACCTAAGACAATCCTGCCTGAAACAGGGCCTGGCATCCCATCAGTGCTCAATAAACA     FFFFF:F:FFFFFFFFF:FFFF,,FF:FFFFFFFFFF:FF:FFF,FFFFF:FFFFF,FFF,FFF,FF,,FF,FFFFFFFFFFFF
#NM:i:0  MD:Z:84 AS:i:84 XS:i:84 XA:Z:14,-19785603,84M,0;       CR:Z:AAGCCGCGTTGACGGA    CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:ACCCTCCT   QT:Z:,F::FFF:   GP:i:698302392 MP:i:698302766   RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1

#A00228:250:H3TVLDSXX:4:1315:24505:4507  1123    22      16057320        27      84M     =       16057594        374    GTAAAATGTAGATTTTGGCACACGGGGACCTAAGACAATCCTGCCTGAAACAGGGCCTGGCATCCCATCAGTGCTCAATAAACA     FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
#NM:i:0  MD:Z:84 AS:i:84 XS:i:84 XA:Z:14,-19785603,84M,0;       CR:Z:AAGCCGCGTTGACGGA    CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:AAGCCGCGTTGACGGA-1 BC:Z:GTTGCAGC   QT:Z:FFFFFFFF   GP:i:698302392 MP:i:698302766   RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:4


# a different set: here, a position duplicate from another CELL is NOT marked duplicate
#A00228:250:H3TVLDSXX:1:1146:8477:8469   1171    22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    NM:i:1  MD:Z:34G65      AS:i:95XS:i:90  XA:Z:14,+19789081,100M,2;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFFFFFF   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:1:1146:8160:8484   147     22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FFFFFFF,FFFFFFFF:FFFFFF:FF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF    NM:i:1  MD:Z:34G65      AS:i:95XS:i:90  XA:Z:14,+19789081,100M,2;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFFFFFFFFFFFFFF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFFFFFF   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:1:2628:1479:16297  1171    22      16053829        27      100M    =       16053921        -8     TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCAAGGACCATGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     FF,F:FFFF:FFFF:FF,,FF,FFFFFFFF,FFFFFFFFFFF,:::FFFFFFFF,,FF::FFF:FFFFFFFFFFFFFFFFFFFFFFFFF,,FFFF:FFFF    NM:i:2  MD:Z:34G7G57    AS:i:90XS:i:85  XA:Z:14,+19789081,100M,3;       CR:Z:ACAGCTAAGGCATTGG   CY:Z:FFFF,FFFFF,FF:FF
#CB:Z:ACAGCTAAGGCATTGG-1 BC:Z:ACCCTCCT   QT:Z:FFFF,F,F   GP:i:698299001  MP:i:698298993  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:1
#
#A00228:250:H3TVLDSXX:2:1373:28854:34428 147     22      16053829        50      100M    =       16053669        -260   TACCAGAGGCATGGGGTGAGGCATGTTGCAGGTCGAGGACCAGGGCCATCTCACTGCCTGAGCCCATGGACTGGCTCAGGGGTCTGTCAGATGATTCTAG     :FFF:FFFFFF::FFF:FFF,FFF:FFFFFFFFFF:FFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFF    NM:i:0  MD:Z:100        AS:i:10XS:i:95  XA:Z:14,+19789081,100M,1;       CR:Z:CGGAGCTAGATCCGAG   CY:Z:FFF:FFFFFFFFFFFF
#CB:Z:CGGAGCTAGATCCGAG-1 BC:Z:CAATGGAG   QT:Z::::FF,,,   GP:i:698299001  MP:i:698298741  RG:Z:bj_mkn45_10pct:MissingLibrary:1:H3TVLDSXX:2










                ## examine every fragment that might contribute to bin
                #foreach my $x(@frags){
                #    
                #    # these fragments are entirely within this bin, before the break
                #    if($$x[1] < $brkStart){ 
                #        $binCount{$ALL}++; # just count them
                #        $binCount{$$x[0]}++;
                #        
                #    # these fragments are likely to be split by the break    
                #    } else { 
                #        my $inc = 1 / $$x[8];
                #        my $CB = $$x[0];
                #        my ($minPos, $maxPos);
                #        # fragment must have come BEFORE the index fragment 
                #        if($$x[3] < $brkStart){ 
                #            my $bin = ($brkStart - $$x[3]) * $inc; 
                #            $binCount{$ALL} += $bin;
                #            $binCount{$CB} += $bin; # these bases always in this bin
                #            $brkTotal{$ALL} += 1 - $bin;
                #            $brkTotal{$CB} += 1 - $bin;
                #            ($minPos, $maxPos) = ($brkStart, $$x[1]); # these bases in break region
                #        
                #        # fragment must have come AFTER the index fragment
                #        } elsif($$x[1] > $brkEnd){ 
                #            my $next = ($$x[1] - $brkEnd) * $inc;
                #            $nextBinCount{$ALL} += $next; # these bases always in next bin
                #            $nextBinCount{$CB} += $next;
                #            $brkTotal{$ALL} += 1 - $next;
                #            $brkTotal{$CB} += 1 - $next;
                #            ($minPos, $maxPos) = ($$x[3], $brkEnd); # these bases in break region
                #
                #        # some fragments may reside entirely in the break region
                #        } else {
                #            $brkTotal{$ALL}++;
                #            $brkTotal{$CB}++;
                #            ($minPos, $maxPos) = ($$x[3], $$x[1]);
                #        }
                #        
                #        # build a coverge map of the break region
                #        # TODO: make this go in 10 bp increments? would reduce loop time
                #        foreach my $pos($minPos..$maxPos){ 
                #            $pos{$ALL}[$pos-$brkStart] += $inc;
                #            $pos{$CB}[$pos-$brkStart] += $inc;
                #        }                           
                #    }
                #}
                #
                ## step through the break region until the definitive breaking base is found
                #foreach my $i(0..$#{$pos{$ALL}}){
                #    foreach my $CB(@allCells){
                #        $brkCumSum{$CB} += $pos{$CB}[$i] || 0;
                #    }
                #    if($binCount{$ALL} + $brkCumSum{$ALL} >= $binWeight){
                #        $binEnd = $brkStart + $i;
                #        last;  
                #    }
                #}