
# get passed agruments
args <- as.list(Sys.getenv())

# load common functions and parameters
command <- "merge10X"
source(paste(args$LIB_DIR, "common.R", sep="/"))
source(paste(args$LIB_DIR, "common10X.R", sep="/"))

# load bins x cells as data frame
read.bgz("bins", args$SAMPLE, header=TRUE)
allI <- 5 # file format is chrom, start, end, nGapBases, all_cells, cell1, cell2, ...

# calculate bin metrics
message("calculating bin metrics")
bin         <- list()
bin$N       <- nrow(d)
bin$sizes   <- (d$end - d$start) - d$nGapBases # the number of ACGT(!N) bases that contributed to the bin
bin$medianSize <- median(bin$sizes)
numChroms  <- sub("CHR", "", toupper(d[[1]])) # make sure chroms are 1,2...X,Y
autosomeIs <- suppressWarnings(!is.na(as.numeric(numChroms)))
bin$medianSizeAutosome <- median(bin$sizes[autosomeIs])
bin$data    <- d[,1:allI]

# calculate cell metrics
message("calculating cell metrics")
cell <- list()
cell$allIs  <- allI:ncol(d)
cell$dIs    <- (allI+1):ncol(d)
cell$Is     <- cell$dIs - allI
cell$N      <- ncol(d) - allI
cell$allSum <- sum(d[[allI]]) # total number of reads over all cells
cell$sums   <- colSums(d[,cell$dIs])  # total number of reads per cell
cell$avgSum <- cell$allSum / cell$N

# report various values
reportCount(bin$medianSize, "median bin size")
reportCount(bin$medianSizeAutosome, "median autosome bin size")
reportCount(cell$allSum, "read pairs over all cells")
reportCount(cell$avgSum, paste("read pairs per each of", commifyInt(cell$N), "normalized cells", sep=" "))
reportCount(cell$avgSum/bin$N, paste("read pairs per each of", commifyInt(bin$N), "bins per cell", sep=" "))

# create sample data array: row=bin(x|i), column=cell(y|i), layer(z)=data
message("assembling data array")
layer <- list()
layer$names <- list(raw="raw", # actual read count for each cell-bin
                     cn="cn",  # adjusted estimated copy number for each cell-bin
                    cnc="cnc", # estimated copy number change relative to the modal copy number 
                    exp="exp") # expected read count to match the modal copy number
layer$N <- length(layer$names)
layers <- array(0, dim = c(bin$N, cell$N, layer$N),
              dimnames = list(bin=NULL, cell=NULL, layer=layer$names))
layers[,,layer$names$raw] <- as.matrix(d[,cell$dIs])

# examine the count range for every cell
message("scaling each cell to quantal read numbers")
minBinN  <- 200 # merge bins until they have at least this many average counts in a cell
searchCN <- 1:20 # should not need to be this large when not using cancer cells with high copy amplifications (or, just mask bins?)
optimizeReadN <- function(readsPerAllele){ # sum of a cell's bin distances to their nearest integer copy number
    sum(sapply(collapsed, function(N) min(abs(N / readsPerAllele - searchCN))))
}
cell$readsPerAllele <- sapply(cell$Is, function(j){ # how many reads in a bin corresponds to one allele for a cell
    #if(j %% 100 == 1 & j != 1) message(paste(j-1, "of", commifyInt(cell$N), "cells completed"))
    binAvg <- cell$sums[j] / bin$N  # TODO: this should probably only use autosomes? #################
    nCollapseBins <- ceiling(minBinN / binAvg) # merge adjacent bins as needed to get the cell's count up 
    collapsed <<- collapseVector(layers[,j,layer$names$raw], nCollapseBins)
    readsPerAlleleGuess <- binAvg * nCollapseBins / modalCNGuess # currently based on default number of alleles
    optimize(optimizeReadN, c(0.5, 1.5) * readsPerAlleleGuess)$minimum / nCollapseBins
})
layers[,,layer$names$cn] <- t(apply(layers[,,layer$names$raw], 1, '/', cell$readsPerAllele)) # STEP 1 of CN calc

# examine the count range for every bin
message("adjusting bins for outlier cells") # this does not mean "removing bad/noisy cells"! 
n20   <- round(cell$N / 5, 0)               # it means excluding possible CNV cells when adjusting to bin CN
mid60 <- n20:(cell$N - n20) # exclude the top and bottom 20% of cell values when optimizing bins
optimizeBinAdj <- function(binAdj){ # sum of middle 60% cells distance to CN=2
    sum(abs(midCN * binAdj - 2))
}
bin$adjustments <- sapply(1:bin$N, function(i){
    #if(i %% 1000 == 1 & i != 1) message(paste(i, "of", commifyInt(bin$N), "bins completed"))
    midCN <<- sort(layers[i,,layer$names$cn])[mid60]
    optimize(optimizeBinAdj, c(1/3, 3))$minimum    
})
layers[,,layer$names$cn] <- apply(layers[,,layer$names$cn], 2, '*', bin$adjustments) # STEP 2 of CN calc

# determine modal copy numbers from bin sizes and outlier adjustments
message("adjusting bins for modal copy number across all cells")
bin$modalCNs <- bin$medianSizeAutosome / bin$sizes * modalCNGuess / bin$adjustments # larger bins have lower copy number across all samples
# TODO: this is where the 50% midpoint is enforced (via round); in fact poor coverage CN2 bins extend a bit past the midpoint
bin$modalCNsInt <- pmax(1, round(bin$modalCNs, 0)) # most likely integer copy number of the predominant cell state
layers[,,layer$names$cn] <- apply(layers[,,layer$names$cn], 2, '*', bin$modalCNsInt / modalCNGuess) # STEP 3 of CN calc
#writeBgz(cn, CN_FILE)

# calculate copy number changes per bin per cell
message("calculating copy number changes") # would this be better as bin$modalCNs?? (probably not)
layers[,,layer$names$cnc] <- apply(layers[,,layer$names$cn], 2, '-', bin$modalCNsInt) 
#writeBgz(cnc, CNC_FILE)

# calculate read counts expected per bin to achieve the modal copy number (used in HMM)
message("calculating expected read counts per cell per bin") 
layers[,,layer$names$exp] <- sapply(cell$readsPerAllele, '*', modalCNGuess / bin$adjustments)

# save object data for HMM
file <- get.outFile("data", args$SAMPLE, "RData")
save.image(file)

# plot value histograms
message("plotting cell read count histogram")
plotHistogram(cell$sums, "readCounts", "Cell Total Read Count", breaks=25)

message("plotting cell reads per allele histogram")
plotHistogram(cell$readsPerAllele, "readPerAllele", "Cell Allele Read Count", breaks=25)

message("plotting bin size histogram")
plotHistogram(bin$sizes, "binSizes", "Bin Size (bp)",
              c(0, 3 * bin$medianSize), bin$medianSizeAutosome)

message("plotting bin adjustment histogram")
plotHistogram(bin$adjustments, "binAdjustments", "Bin Adjustment")

message("plotting bin modal copy number histogram")
plotHistogram(bin$modalCNs, "modalCopyNumbers", "Bin Modal Copy Number",
              c(0, 5), 2)

message("plotting bin copy number histogram across all cells")
plotHistogram(as.vector(layers[,,layer$names$cn]), "observedCopyNumbers", "Bin Copy Number",
              c(0, 5), 2)

message("plotting copy number change histogram")
plotHistogram(as.vector(layers[,,layer$names$cnc]), "copyNumberChange", "Bin Copy Number Change",
              c(-3, 3), 0)

#
## make CNC heat maps
#message("plotting chromosome heat maps")
#cnColors <- c("red","orange","yellow")
##pd <- t(layers[,,layer$names$cnc])
##pd <- ifelse(pd >  2,  2, pd) # make symmetric and center for heat map
##pd <- ifelse(pd < -2, -2, pd)
##png(filename = paste(PLOT_DIR, "/", SAMPLE, ".GENOME.cnc.png", sep=""),
##    width = plotDim, height = plotDim, units = "px", pointsize = pointsize,
##    type = c("cairo"))
##hmp <- heatmap(pd, Colv=NA, scale="none", ColSideColors=cnColors[modalCNsInt],
##        main="GENOME", xlab="Bins", ylab="Cells",
##        labRow=rep("",nrow(pd)), labCol=rep("",ncol(pd)) )    
##graphics.off()
#for(chrom in unique(d[[1]])){
#    message(paste("   ", chrom))
#    binIs <- which(d[[1]]==chrom)
#    pd <- t(layers[binIs,,layer$names$cnc])
#    pd <- ifelse(pd >  2,  2, pd) # make symmetric and center for heat map
#    pd <- ifelse(pd < -2, -2, pd)
#    file <- get.plotFile("heatmap", paste(args$SAMPLE, chrom, sep="."))
#    png(file, width=plotDim, height=plotDim,
#        units="px", pointsize=pointsize, type = c("cairo")) 
#    heatmap(pd, Colv=NA, scale="none", ColSideColors=cnColors[bin$modalCNsInt[binIs]],
#            main=chrom, xlab="Bins", ylab="Cells",
#            labRow=rep("",nrow(pd)), labCol=rep("",ncol(pd)) )
#    graphics.off()
#}
