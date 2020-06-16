#!/bin/bash

#q    require $SAMPLE $MASTERS_DIR $COV_DIR $SEG_DIR 
#q    require $BIN_SIZE 

#$    -N  optimize_$SAMPLE
#$    -wd $SEG_DIR
#$    -l  vf=4G

#sleep 30 

SLP_GLOB="$SEG_DIR/$SAMPLE.slope_*.$BIN_SIZE.bed.bgz"
COV_BED="$COV_DIR/$SAMPLE.binCoverage.$BIN_SIZE.bed.bgz"

echo "optimizing sloped segments for $SAMPLE bin=$BIN_SIZE MIN_BIN_READS=$MIN_BIN_READS"

gunzip -c $SLP_GLOB |
cut -f1-8 |
awk '{print $0"\t"$2}' | 
bedutil split |
sort -k1,1 -k2,2n |
bedtools intersect -loj -a - -b $COV_BED |
cut -f 1-9,11,14,16,17 |
groupBy -g 1,2,3,4,5,6,7,8,9 -c 10,11,12,13 -o collapse,collapse,collapse,collapse |
groupBy -g 1,2,3 -c 7,8,9,10,11,12,13 -o collapse,collapse,collapse,distinct,distinct,distinct,distinct | 



head -n 100
exit 1


# collect the bin coverage data from this chromosome
slurp tabix $COV_BED $CHROM |
# smooth the RFD trace, i.e. remove high frequency noise
# disregard bins with too low or too high of a total read count to be reliable
awk 'BEGIN{OFS="\t"}($7+$8)>'$MIN_BIN_READS'&&($7+$8)<='$MAX_BIN_READS'{print $5, $0}' |
smooth -j $SMOOTH_J |
# adjust the bin crick and watson counts to match the smoothed trace
# maintain original RFD and read counts for later determination of weighted segment slopes
awk 'BEGIN{
    OFS = "\t";
}{
    if($1 > 1){ $1 = 1} else if ($1 < -1) { $1 = -1 }
    N = $8 + $9;
    c = int( (($1 + 1) * N) / 2 + 0.5 );
    w = N - c;
    print $2,$3,$4,$5,$1,$7,c,w,$6,$8,$9;
}' |
bgzip -c | 
slurp -o $SMOOTHED_BED
checkPipe

sleep 30
echo $SMOOTHED_BED
gunzip -c $SMOOTHED_BED | wc -l

# CODE BELOW IS THE CURRENT SEGMENT THRESHOLD APPROACH

echo "parsing segments at threshold=$SEGMENT_THRESHOLD"

# collect all bin coverage data, merge to smoothed RFD and mark gap bins
slurp tabix $COV_BED $CHROM | # :11620323-12043105
bedtools intersect -loj -a stdin -b $SMOOTHED_BED |
cut -f 1-9,13 |
bedtools intersect -c -a stdin -b $GAP_FILE |
# find ascending and descending segments
perl $MASTERS_DIR/slaves/parse_segments.pl $SEGMENT_THRESHOLD $FUSE_MAX_FRAC $FUSE_MAX_BINS |
# calculate linear regression of RFD over each called segment
groupBy -g 1 -c 4,7,9,10,2,4,5 -o collapse,collapse,collapse,collapse,distinct,min,max |
Rscript $MASTERS_DIR/slaves/slope_slopes.R |
# commit the final results to BED file
# format: chrom, start, end, segmentI, sloped(0=~flat), strand=sign of slope, slope, intercept, pos, resd
awk 'BEGIN{OFS="\t"}{print "'$CHROM'", $7, $8, $1, $6, ($3 > 0) ? "+" : "-", $3, $2, $4, $5}' |
bgzip -c | 
slurp -o $SLP_BED
checkPipe

waitForFile $SLP_BED
echo $SLP_BED
gunzip -c $SLP_BED | wc -l

echo "done"


