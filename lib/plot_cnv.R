
# get passed arguments
plotFile     <- Sys.getenv("plotFile")
spanFile     <- Sys.getenv("spanFile")
snpFile      <- Sys.getenv("snpFile")
gapDupFile   <- Sys.getenv("gapDupFile")
plotDir      <- Sys.getenv("plotDir")
centerWidth  <- as.numeric(Sys.getenv("centerWidth"))
maxSnpCount  <- as.numeric(Sys.getenv("maxSnpCount"))
refName      <- Sys.getenv("refName")
altName      <- Sys.getenv("altName")

# set variables
widthPadding <- (centerWidth - 1) / 2;
idx_col <- 'darkred'
oth_col <- 'blue'
non_col <- 'grey'
cnlim   <- c(-2,4)
snplim  <- c(-maxSnpCount,maxSnpCount)

# get the data
write("  loading data", file=stderr())
d <- read.table(plotFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(d) <- c('chrom', 'region', 'group', 'grpId', 'cnvType', 'cnvSize', 'cnvFrac',
                 'idx_sample', 'oth_sample',
                 'start', 'end', 'left', 'right',
                 'sample', 'center', 'modCN', 'dCN');
s <- read.table(spanFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(s) <- c('grpId', 'idx_sample', 'sample', 'cnvType', 'start', 'end', 'row');
S <- read.table(snpFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(S) <- c('grpId', 'idx_sample', 'sample', 'pos', 'nRef', 'nAlt', 'refBase', 'altBase');
g <- read.table(gapDupFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(g) <- c('grpId', 'start', 'end', 'col');

# plot functions
addVLines <- function(){
    abline(v=c(left, right)/1e6, col="grey", lty=2)
}
addGaps <- function(lim){
    if(nrow(gg) > 0){
        for(i in 1:nrow(gg)){
            rect(gg[i,'start']/1e6, lim[1], gg[i,'end']/1e6, lim[2], col=gg[i,'col'], border=NA)
        }
    }
}
plotCopyNumber <- function(){
    legend("topleft", c('Modal Copy Number', idx_sample), col=c('black', idx_col), lwd=1, bty="n", cex=1)
    addGaps(cnlim)
    abline(h=seq(-2,4,1), col="grey", lty=2)    
    abline(h=0, col="black", lty=1)
    addVLines()
    isFirst <- TRUE
    for (sample in unique(ddd[,'sample'])){    
        dddd <- ddd[ddd$sample==sample,]
        if(nrow(dddd)>0){
            if(isFirst){
                dddd[,'modCN'] <- ifelse(dddd$modCN < cnlim[1], cnlim[1], dddd$modCN)
                dddd[,'modCN'] <- ifelse(dddd$modCN > cnlim[2], cnlim[2], dddd$modCN) 
                lines (dddd$center/1e6, dddd$modCN, col="black")
                isFirst <- FALSE
            }                 
            dddd[,'dCN'] <- ifelse(dddd$dCN < cnlim[1], cnlim[1], dddd$dCN)        
            dddd[,'dCN'] <- ifelse(dddd$dCN > cnlim[2], cnlim[2], dddd$dCN)
            if(sample==idx_sample) {   
                lines (dddd$center/1e6, dddd$dCN, col=idx_col)
                points(dddd$center/1e6, dddd$dCN, col=idx_col, pch=20)
            } else if(sample==oth_sample) {
                lines (dddd$center/1e6, dddd$dCN, col=oth_col)   
            } else {
                lines (dddd$center/1e6, dddd$dCN, col=non_col)
            }
        }
    }  
}
getBias <- function(snps){
    nRef <- sum(snps$nRef)
    nAlt <- sum(snps$nAlt)
    nReads <- nRef + nAlt
    ifelse(nReads>0, (nRef - nAlt) / nReads * snplim[2], 0)
}
plotSnpCounts <- function(){
    addGaps(snplim)
    abline(h=((2 - 1) / 3 * snplim[2]), col="grey", lty=1)
    abline(h=((1 - 2) / 3 * snplim[2]), col="grey", lty=1)
    abline(h=((1 - 0) / 1 * snplim[2]), col="grey", lty=1)
    abline(h=((0 - 1) / 1 * snplim[2]), col="grey", lty=1)
    abline(h=0, col="black", lty=1)
    addVLines()
    mtext(refName, side=2, line=2, at=snplim[2]*0.9, las=1, cex=0.8)
    mtext(altName, side=2, line=2, at=snplim[1]*0.9, las=1, cex=0.8)
    SSS <- SS[SS$sample==idx_sample,]    
    if(nrow(SSS)>0){
        nRef <- ifelse(SSS$nRef > snplim[2],  snplim[2],  SSS$nRef)                 
        nAlt <- ifelse(SSS$nAlt > snplim[2], -snplim[2], -SSS$nAlt)        
        points(SSS$pos/1e6, nRef, col="black", type="h")
        points(SSS$pos/1e6, nAlt, col="black", type="h")        
    }
    for(sample in c(oth_sample, idx_sample)){
        SSS <- SS[SS$sample==sample,]
        if(nrow(SSS)>0){
            cnvBias <- getBias(SSS[SSS$pos>=left&SSS$pos<=right,])
            lftBias <- getBias(SSS[SSS$pos>=start&SSS$pos<left,])
            rgtBias <- getBias(SSS[SSS$pos>right&SSS$pos<=end,])
            col <- if(sample==oth_sample) {oth_col} else {idx_col}
            lines(c(left,  right) / 1e6, rep(cnvBias,2), col=col, lwd=2)
            lines(c(start, left)  / 1e6, rep(lftBias,2), col=col, lwd=2)
            lines(c(right, end)   / 1e6, rep(rgtBias,2), col=col, lwd=2)                 
        }
    }   
}
plotSpans <- function(sample, col){
    legend("topleft", sample, col=col, lwd=1, bty="n", cex=1)
    addGaps(ylim)
    addVLines()    
    sss <- ss[ss$sample==sample,]
    if(nrow(sss)>0){
        segments(sss$start/1e6, sss$row, sss$end/1e6, sss$row, col=col, lwd=1, lend=2)
        abline(h=0.5, col="grey", lty=2)
        for(cnvType in c('del', 'dup', 'invF', 'invR')){
            ssss <- sss[sss$cnvType==cnvType,]
            if(nrow(ssss)>0){
                min <- min(ssss[,'row'])
                max <- max(ssss[,'row'])
                abline(h=max+0.5, col="grey", lty=2)
                mtext(cnvType, side=2, line=1, at=min+(max-min)/2, las=1, cex=0.8)                
            }    
        }   
    }
}

# make the plots
for (grpId in unique(c(d$grpId, s$grpId))){
    
    # get this index event
    write(paste("  ", grpId, sep=""), file=stderr())    
    dd      <- d[d$grpId==grpId,]
    ss      <- s[s$grpId==grpId,]
    SS      <- S[S$grpId==grpId,]
    gg      <- g[g$grpId==grpId,]
    group   <- dd[1,'group']
    
    # set the coordinates
    chrom   <- dd[1,'chrom']
    start   <- dd[1,'start']
    end     <- dd[1,'end']
    left    <- dd[1,'left']
    right   <- dd[1,'right']
    clim    <- c(start, end) /  1e6

    # determine the index and other samples
    idx_samples <- unique(c(dd$idx_sample, ss$idx_sample)) # can be more than one if region has multiple candidate samples
    idx_sample  <- idx_samples[1]
    ddd <- dd[dd$idx_sample==idx_sample,] # only need to plot the trace from one call
    oth <- ddd[ddd$sample!=idx_sample&ddd$center>=left&ddd$center<=right,]    
    dlt <- abs(oth[,'dCN'])    
    oth_sample <- if(!is.na(ddd[1,'oth_sample'])){
        ddd[1,'oth_sample']
    } else if(length(dlt)>0){
        oth[which(dlt==max(dlt))[1],'sample']
    } else {
        'xxxx'
    }

    # open the jpg file
    jpgFile <- paste(ddd[1,'region'], group, grpId, "jpg", sep=".")
    jpgFile <- paste(plotDir, jpgFile, sep="/")        
    jpeg(jpgFile,
         width=5, height=6, units="in", pointsize=9,
         res=400, quality=80)    
    
    # set the composite layout
    layout(matrix(c(1,2,3,4), 4, 1, byrow=TRUE), heights=c(0.28,0.22,0.22,0.28))

    # plot the copy number traces
    par(mar=c(0.1, 4.1, 4.1, 0.5))
    plot(0, 0, type="n", xlim=clim, ylim=cnlim,
         xaxt="n",
         ylab="Copy Number (Change)",
         main=paste(paste(group, " #", grpId, ", ",
                          ddd[1,'cnvType'], " ", as.numeric(ddd[1,'cnvSize']), " bp", sep=""))
    )
    plotCopyNumber()

    # plot the SNP counts
    par(mar=c(0.1, 4.1, 0.1, 0.5))    
    plot(0, 0, type="n", xlim=clim, ylim=snplim,
         xaxt="n",
         ylab="SNP Count"
    )
    plotSnpCounts()

    # plot the index sample spans
    par(mar=c(0.1, 4.1, 0.1, 0.5))
    rows <- ss[ss$sample==idx_sample,'row']
    ylim <- c(0, 1 + if(length(rows)>0) {max(rows)} else {0})
    plot(0, 0, type="n", xlim=clim, ylim=ylim,
         xaxt="n",
         yaxt="n", ylab=""
    )
    plotSpans(idx_sample, idx_col)
    
    # plot the other/control sample spans
    par(mar=c(4.1, 4.1, 0.1, 0.5))
    rows <- ss[ss$sample==oth_sample,'row']
    ylim <- c(0, 1 + if(length(rows)>0) {max(rows)} else {0})
    plot(0, 0, type="n", xlim=clim, ylim=ylim,
         xlab=paste(chrom, "coordinate (Mb)", sep=" "),
         yaxt="n", ylab=""
    )
    plotSpans(oth_sample, oth_col)

    # finish the jpg
    graphics.off()
}
