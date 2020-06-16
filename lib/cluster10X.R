
# get passed agruments
clusterArgs <- as.list(Sys.getenv())
args <- clusterArgs

# load common functions and parameters
clusterCommand <- "cluster10X"
command <- clusterCommand
source(paste(clusterArgs$LIB_DIR, "common.R", sep="/"))
source(paste(clusterArgs$LIB_DIR, "common10X.R", sep="/"))

# load the source data as output by svtools merge10X
message(paste("loading sample:", clusterArgs$SAMPLE, sep=" "))
load(get.inFile("data", clusterArgs$SAMPLE, "RData", 'merge10X'))
args <- clusterArgs # since loaded file had its own args... (args object used by common functions)
command <- clusterCommand

# remove the big and unneeded d object for memory efficiency...
#   bin 2.8 Mb <<<< sizes of the biggest loaded objects from one example file
#   d 184.6 Mb <<<< remove this immediately
#   layers 732.6 Mb <<<< remove this prior to saving Z score layers
rm(d)    

# disregard potentially untrustworthy bins
# good bins are those with reliable data counts that don't cross large gaps
# TODO: more aggressive bin quality thresholds? would need mappability most likely
goodBins <- bin$modalCNs > 0.5 & bin$data$nGapBases <= 1000

# disregard potentially untrustworthy cells
# currently based solely on CellRanger "noisy" call
# TODO: more aggressive cell quality thresholds (hyper segmented?)
cr_stats <- read.csv(paste(args$CR_DIR, 'per_cell_summary_metrics.csv', sep="/"),
                      header=TRUE, stringsAsFactors=FALSE)
goodCells <- (1:cell$N)[cr_stats$is_noisy==0]

print(goodCells)
quit('no')

# get output chromosomes; save genome info for server
genome <- clusterArgs$GENOME
chromI <- 1
chroms <- getOrderedChroms(bin$data[[chromI]], as.logical(clusterArgs$INCLUDE_Y))
save(genome, chroms, file = get.outFile("chroms", clusterArgs$SAMPLE, "RData"))

# set some processing and output variables
zLN <- list( # names of the output layers in Z score array
  cn="cn",
  raw="raw",
  exp0="exp0",
  expG="expG",
  expL="expL",
  z0="z0",
  zG="zG",
  zL="zL",
  HMM="HMM"
)
minExp <- 0.01
collapseFactor <- 10 # for whole genome view
HMM_set_persistence(1 - 1e-6)

# create a layers array with bin counts and Z-scores
assembleArray <- function(bins, collapseFactor=NULL){
  
  # output only has good bins, must account for this in future coordinate-based plots
  binIs <- which(bins & goodBins)
  
  # set base counts for each bin x cell combination
  # collapse (i.e. merge) adjacent bins as requested
  if(is.null(collapseFactor)){
    zL <- array(0, dim = c(length(binIs), cell$N, length(zLN)),
                dimnames = list(bin=NULL, cell=NULL, layer=zLN))
    zL[,,zLN$cn]   <- layers[binIs,,layer$names$cn]
    zL[,,zLN$raw]  <- layers[binIs,,layer$names$raw]
    zL[,,zLN$exp0] <- layers[binIs,,layer$names$exp]
    inc <- apply(zL[,,zLN$exp0], 2, '*', 1 / bin$modalCNsInt[binIs])
  } else {
    mCN <- round(collapseVector(bin$modalCNsInt[binIs], collapseFactor) / collapseFactor, 0)
    zL  <- array(0, dim = c(length(mCN), cell$N, length(zLN)),
                 dimnames = list(bin=NULL, cell=NULL, layer=zLN))
    zL[,,zLN$cn]   <- apply(layers[binIs,,layer$names$cn],  2, collapseVector, collapseFactor) / collapseFactor
    zL[,,zLN$raw]  <- apply(layers[binIs,,layer$names$raw], 2, collapseVector, collapseFactor)
    zL[,,zLN$exp0] <- apply(layers[binIs,,layer$names$exp], 2, collapseVector, collapseFactor)
    inc <- apply(zL[,,zLN$exp0], 2, '*', 1 / mCN)
  }
  
  # calculate CNV counts and Z scores for each bin x cell combination
  zL[,,zLN$expG] <- zL[,,zLN$exp0] + inc
  zL[,,zLN$expL] <- zL[,,zLN$exp0] - inc
  zL[,,zLN$expL] <- ifelse(zL[,,zLN$expL] < minExp, minExp, zL[,,zLN$expL]) # prevent divide by zero
  zL[,,zLN$z0]   <- (zL[,,zLN$raw] - zL[,,zLN$exp0]) / sqrt(zL[,,zLN$exp0])
  zL[,,zLN$zG]   <- (zL[,,zLN$raw] - zL[,,zLN$expG]) / sqrt(zL[,,zLN$expG])
  zL[,,zLN$zL]   <- (zL[,,zLN$raw] - zL[,,zLN$expL]) / sqrt(zL[,,zLN$expL])
  
  # return the array
  zL
}

# print an RData object for segmentation and R Shiny visualization
saveZLayers <- function(chrom){
  name <- paste(clusterArgs$SAMPLE, chrom, sep=".")
  objs <- ls(env=globalenv())
  objs <- objs[!(objs %in% "layers")] # remove prior layers object for memory
  save(list = objs,
       file = get.outFile("Z_scores_HMM", name, "RData"))   
}

# hierachically cluster the cells
message("clustering cells over entire genome (collapsed bins)")
zLayers <- assembleArray(rep(TRUE, bin$N), collapseFactor) # will be a new array for each chrom
hClustGenome <- hclust(pearson.dist(t(zLayers[,,zLN$z0])))   # will save as is for all chromosomes also
zLayers[,,zLN$HMM] <- sapply(1:cell$N, HMM_viterbi)
saveZLayers('all')

# cluster each chromosome
# segment individual cells (will work on clustered cell HMM later)
message("clustering cells over each chromosome individually (uncollapsed bins)")
for (chrom in chroms){
  message(paste(" ", chrom))
  chromBins <- bin$data[[chromI]] == chrom
  zLayers <- assembleArray(chromBins) 
  hClustChrom <- hclust(pearson.dist(t(zLayers[,,zLN$z0])))
  # hc <- hclust(pearson.dist(t(rbind(cL[,,zLN$z0],cL[,,zLN$zG],cL[,,zLN$zL])))) # or just z0?
  zLayers[,,zLN$HMM] <- sapply(1:cell$N, HMM_viterbi)
  saveZLayers(chrom)
}

