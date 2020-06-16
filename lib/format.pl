use strict;
use warnings;

#========================================================================
# 'format.pl' has subs for handling file headers and columns
#========================================================================

#========================================================================
# define variables
#------------------------------------------------------------------------
use vars qw($utility $command $slurp $error);
#========================================================================

#========================================================================
# read headers from 'extract' or 'map' output files (gz, per sample+library)
#------------------------------------------------------------------------
sub getFileHeader {
    my ($strpHdr, $type, @required) = @_;
	my %tags = map { if($_ =~ m/(.+):(.+)/){
						 $1 => $2
					 } else {
						 "x" => 1
					 } } split("\t", $strpHdr);
	($tags{application} and $tags{application} eq "svtools" and
     $tags{file_type}   and $tags{file_type}   eq $type)
        or die "$error: invalid $type header: #$strpHdr";
    foreach my $tag(@required){
        $tags{$tag} or die "$error: could not find '$tag' in $type header: #$strpHdr";
    }
	return %tags;
}
#========================================================================

#========================================================================
# 'find' file cols
#------------------------------------------------------------------------
my $findI = 0;
our %findCols = map { $_ => $findI++ } qw(
  chrom1 minPThis maxPThis svID svSze strand
  chrom2 minPThat maxPThat
  jxnSource jxnType
  frags
  nPosSmp nCllSmp
  fracOvlp fracOvlp1 fracOvlp2
  pSet pOutlier    
); # samples evidence (all, s1, s2 ...)
our $nFindLeadCol = $findCols{pOutlier} + 1;
#------------------------------------------------------------------------
my $cndFindI = 0;
our %cndFindCols = map { $_ => $cndFindI++ } qw(
  chrom start end svID svSze strand
  jxnSource jxnType
  frags
  nPosSmp nCllSmp
  fracOvlp fracOvlp1 fracOvlp2
  pSet pOutlier    
); # samples evidence (all, s1, s2 ...)
our $nCndFindLeadCol = $cndFindCols{pOutlier} + 1;
#========================================================================

#========================================================================
# tabix indexed file cols
#------------------------------------------------------------------------
my $indxI = 0;
our %indxCols = map { $_ => $indxI++ } qw(
  chrom1 minPThis maxPThis svID svSze strand
  chrom2 minPThat maxPThat
  jxnSource jxnType
  nPosSmp nCllSmp
  fracOvlp fracOvlp2 fracOvlp2
  pSet pOutlier    
);
#========================================================================

#========================================================================
# 'call' file cols
#------------------------------------------------------------------------
our %binCols = ( # as assembled in memory, i.e. no chrom
    chrom       => -1,      
    start       =>  0,
    end         =>  1, 
    name        =>  2, # unused
    length      =>  3,
    strand      =>  4, # unused
    noGapLen    =>  5, 
    medCov      =>  6,
    adjLen      =>  7,
    cnRaw       =>  8,
    cnHMM       =>  9,
    cnvCalls    =>  10, # number of samples with any CNV call
    # samples raw counts
    # samples normalized counts
    # samples delta copy number raw
    # samples delta copy number HMM
    # samples cnv call (-1 or 1 for loss/gain)
);
our %binColOffsets = (
    smpCovRaw   =>  0,
    smpCovNrm   =>  1,
    smpDCNRaw   =>  2,
    smpDCNHMM   =>  3,
    smpCNVCall  =>  4,
);
our $nBinLeadCol = $binCols{cnvCalls} + 1;
our %binFileCols = map { $_ => ($binCols{$_} + 1) } keys %binCols;
our $nbinFileLeadCol = $nBinLeadCol + 1;
#------------------------------------------------------------------------
my $cnvI = 0;
our %cnvCols = map { $_ => $cnvI++ } qw(
  chrom start end cnvId nBins strand
  cnvType cnvSize sample sampleIndex
  uniqueness maxNSamples binsList
); # samples x mean(delta_copy_number_raw)
   # samples x PRatio(cnvType vs. neutral)
our $nCnvLeadCol = $cnvCols{binsList} + 1;
#========================================================================

#========================================================================
# 'scan' file cols
#------------------------------------------------------------------------
my $scnPosI = 0;
our %scnPosCols = map { $_ => $scnPosI++ } qw(
  chrom start end sample nFrags strand
  fracLowQual tLens isGap n_log_p deltaCDF_abs deltaCDF_sgn
);
#------------------------------------------------------------------------
my $scnCnvI = 0;
my @scnCnvCols = qw(
  chrom start end sample nBins strand
  nFrags fracLowQual tLens isGap n_log_p deltaCDF_abs deltaCDF_sgn
  deltaCDF_alt cnvType cnvSize cnvFrac
);
our %scnCnvCols = map { $_ => $scnCnvI++ } @scnCnvCols;	
#------------------------------------------------------------------------
my $scnRegI = 0;
our %scnRegCols = map { $_ => $scnRegI++ } qw(
    chrom start end regId nRegBins strand 
	nPosSmp nCllSmp
	idxSmp maxCnvDltIdx maxPosDltIdx maxCnvDltOther maxPosDltOther
	isGap nOtherGap
    idxBin nBins nFrags
	fracLowQual tLens isGap n_log_p deltaCDF_abs deltaCDF_sgn
	deltaCDF_alt cnvType cnvSize cnvFrac
);	
#========================================================================

#========================================================================
# 'snp' file cols
#------------------------------------------------------------------------
my $gntPosI = 0;
our %gntCols = map { $_ => $gntPosI++ } qw(
  chrom start end snpId nReads strand
  refBase altBase nRef nAlt readBases
);
#------------------------------------------------------------------------
my $poolPosI = 0;
our %poolCols = map { $_ => $poolPosI++ } qw(
    chrom start end snpId nReads strand
    refBase altBase nRef nAlt
);
#------------------------------------------------------------------------
my $biasSnpI = 0;
our %biasSnpCols = map { $_ => $biasSnpI++ } qw(
    chrom start end snpId nReads strand
    refBase altBase nRef nAlt biasType pRatio
);
#------------------------------------------------------------------------
my $biasCnvI = 0;
our %biasCnvCols = map { $_ => $biasCnvI++ } qw(
    chrom start end cnvId nSnps strand
    sample nRef nAlt biasType pRatio
);
#------------------------------------------------------------------------
my $biasRegI = 0;
our %biasRegCols = map { $_ => $biasRegI++ } qw(
    chrom start end regId score strand
	nPosSmp nCllSmp
	idxSmp maxPRatioIdx maxPRatioOther
    start_ end_ cnvId nSnps nRef nAlt biasType
	cnvSize bias maxOtherBias
);			
#========================================================================

#========================================================================
# 'plot_prep' standard cols
#------------------------------------------------------------------------
my $pPrpI = 0;
our %pPrpCols = map { $_ => $pPrpI++ } qw(
  chrom start end region score strand
  group grpId cnvType cnvSize cnvFrac sample othSample
);
#========================================================================

1;


##========================================================================
## 'mle' file cols
##------------------------------------------------------------------------
#my $mlePosI = 0;
#our %mlePosCols = map { $_ => $mlePosI++ } qw(
#  chrom start end sample nFrags strand
#  deltaCDF_ref deltaCDF_alt cnvSize cnvFrac cnvType cnvType2
#);
##------------------------------------------------------------------------
#my $mleCnvI = 0;
#my @mleCnvCols = qw(
#  chrom start end sample nBins strand
#  deltaCDF_ref deltaCDF_alt cnvSize cnvFrac cnvType cnvType2
#);
#our %mleCnvCols = map { $_ => $mleCnvI++ } @mleCnvCols;
##------------------------------------------------------------------------
#my $clpRegI = 0;
#our %clpCols = map { $_ => $clpRegI++ } qw(
#    chrom start end cnvId nBins strand 
#	nPosSmp nCllSmp fracGapSmp
#	idxSmp nCnvBins
#    maxDltIdx maxDltOther
#	nFrags deltaCDF_alt cnvSize cnvFrac cnvType cnvType2
#); 
##========================================================================
