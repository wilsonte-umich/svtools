
# get passed agruments
args <- as.list(Sys.getenv())
args$BIN_WEIGHT <- as.numeric(args$BIN_WEIGHT)
args$TRANS_PROB <- as.numeric(args$TRANS_PROB)

# some functions
writeBgz <- function(table, filename){
    message(paste("writing file:", filename, sep=" "))
    bgz <- paste("bgzip -c >", filename, sep=" ")
    write.table(table, file=pipe(bgz), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=FALSE)
    tabix <- paste("tabix -p bed", filename, sep=" ")
    system(tabix)
    print(head(table[,1:10]))
}

## load bins x cells as data frame
#message("loading data")
#message(args$BINS_FILE)
#d <- read.table(args$BINS_FILE, header=TRUE, sep="\t",
#                stringsAsFactors=FALSE, comment.char="",)
#dc <- d[d[[1]]==args$CHROM,]
#print(head(dc[,1:10]))
#message(args$CNC_FILE)
#cnc <- read.table(args$CNC_FILE, header=FALSE, sep="\t",
#                  stringsAsFactors=FALSE, comment.char="",)
#cncc <- cnc[cnc[[1]]==args$CHROM,]
#print(head(cncc[,1:10]))
#rm(cnc)
#
## calculate various needed things
#allN       <- 4 # file format is chrom, start, end, all_cells, cell1, cell2, ...
#allCellIs  <- 4:ncol(d)
#cellIs     <- 5:ncol(d)
#nCells     <- ncol(d) - 4
#nBins       <- nrow(d)
#nBinsC      <- nrow(dc)
#modalCNsInt <- cncc[[allN]]
#allCellSums   <- colSums(d[,allCellIs])  # observed read count for all_cells and each individual cell
#normCellCount <- allCellSums[1] / nCells # total read count assigned to each normalized cell
#rm(d)
#
## normalize the data between samples
#message("normalizing counts between samples")
#dcn <- dc
#for(i in cellIs) dcn[[i]] <- round(dcn[[i]] * normCellCount / allCellSums[i-3], 2)
#normMeds <- apply(dcn[,cellIs], 1, median)
#rm(dcn)
#
## initialize HMM object
#message("loading HMM")
#load(args$HMM_FILE)
#message("setting transition probabilities")
#HMM <- HMM_set_persistence(HMM, 1 - args$TRANS_PROB)


#####################
rdf <- paste(args$BINS_FILE, ".TMP.RData", sep="")
#save.image(rdf)
load(rdf)
HMM <- HMM_set_persistence(HMM, 1 - 1e-6)

# run the HMM on each chromosome
message("solving HMM for each cell individually")
message(paste(nCells, "cells"))
dco <- cncc
for (i in cellIs){
#for (i in cellIs[1]){
    message(paste("  cell", i))
    obs <- HMM_collapse_obs(
        HMM,
        ot=list(
            MEDIAN = round(10 * normMeds * allCellSums[i-3] / normCellCount, 0), # median count of bin a over all cells, scaled to the index cell's total read count
            CN     = modalCNsInt # integer copy number of the bin for most/expected cells
        ),
        os=list(
            N = round(dc[[i]], 0) # integer raw read count for a cell (or combination of cells) in a bin
        )
    )
    dco[[i]] <- CNCs[HMM_viterbi(HMM, obs)$hsi]
}
writeBgz(dco, args$CALL1_FILE)


# general plot settings
plotDim <- 900
pointsize <- 12
cnColors <- c("red","orange","yellow")

pd <- t(as.matrix(cncc[,cellIs]))
pd <- ifelse(pd >  2,  2, pd) # make symmetric and center for heat map
pd <- ifelse(pd < -2, -2, pd)
png(filename = paste(args$PLOT_DIR, "/TESTING.cncc.png", sep=""),
    width = plotDim, height = plotDim, units = "px", pointsize = pointsize,
    type = c("cairo"))
hmr <- heatmap(pd, Colv=NA, scale="none", ColSideColors=cnColors[modalCNsInt],
        main="CNC", xlab="Bins", ylab="Cells", keep.dendro=TRUE,
        labRow=rep("",nrow(pd)), labCol=rep("",ncol(pd)) )    
graphics.off()
    
pd <- t(as.matrix(dco[,cellIs]))
#pd <- ifelse(pd >  2,  2, pd) # make symmetric and center for heat map
#pd <- ifelse(pd < -2, -2, pd)
png(filename = paste(args$PLOT_DIR, "/TESTING.dco.png", sep=""),
    width = plotDim, height = plotDim, units = "px", pointsize = pointsize,
    type = c("cairo"))
heatmap(pd, Rowv=hmr$Rowv, Colv=NA, scale="none", ColSideColors=cnColors[modalCNsInt],
        main="DCO", xlab="Bins", ylab="Cells",
        labRow=rep("",nrow(pd)), labCol=rep("",ncol(pd)) )    
graphics.off()  
    
q("no")

T <- c(rep(0,3), allCellSums / nBins)
for(i in cellIs){
    for(j in i:cellIs){
        Tij <- mean(T[i,j]) # this value changes depending on the pair
        Ni <- dc[,i] / T[i] * Tij # therefore these also change
        Nj <- dc[,j] / T[j] * Tij # therefore the poisson p values will differ based on the read counts in the cells, not ideal


        HMM_dist <- log(pmax(1e-99, dpois(Ns, exp)))
    }
}





#
## update informativity based on model output
## probes in monoallelic runs of homozygosity are never considered informative
## probes in a biallelic state have informativity determined by observed zygosity
#message("calculating derived output columns")
#d$INF_OUT <- mapply(function(IS_NA, AS, INF){
#    if(IS_NA==1){
#        NA
#    } else if(astasms[[AS]]<=1){
#        'no'
#    } else {
#        INF
#    }
#}, d$IS_NA, d$AS_OUT, d$INF_OUT)
#
## calculate changes in output relative to input modeled states
#d$CNC <- if(modeling) NA else ifelse(d$IS_NA, NA, d$CN_OUT - d$CN_IN)
#d$LOH <- if(modeling) NA else ifelse(d$IS_NA, NA, ifelse(d$INF_OUT != d$INF_IN, 'yes', 'no'))
#
## fit and plot refined statistics across array
#cn_col <- 'CN_OUT'
#source(paste(LIB_DIR, "model_array.R", sep="/"))
#
## write the final output probes BED file
#message("writing probes BED")
#header <- paste(PROBE_COL, collapse="\\t")
#header <- paste("#", header, sep="") # comment the header for bgzip
#pipe <- paste("awk 'BEGIN{OFS=\"\\t\";print \"", header,"\"}{$2+=0;$3+=0;print $0}' | bgzip -c > ",
#              PROBES_FILE, sep="")
#write.table(
#    d[,PROBE_COL],
#    file=pipe(pipe),
#    quote=FALSE,
#    sep="\t",
#    row.names=FALSE,
#    col.names=FALSE
#)
#nPrbOut <- nrow(d)
#message(paste("   ", nPrbOut, "probes written to disk"))
#
## write the final output structural variants file as needed
#if(!modeling) source(paste(LIB_DIR, "find_SVs.R", sep="/"))
#
#
##########################
##PROBES_FILE <- '/home/wilsonte_lab/club_house/data/Glover/Irene/samples/msvtools.segment.probes.Scr1-1.bed.bgz'
##RD_FILE <- paste(PROBES_FILE, "RData", sep=".")
###save.image(file=RD_FILE)
###load(RD_FILE)
##source(paste(LIB_DIR, "run_HMM.R", sep="/"))
#



 