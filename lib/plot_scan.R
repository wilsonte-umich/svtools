
# get passed arguments
plotFile     <- Sys.getenv("plotFile")
spanFile     <- Sys.getenv("spanFile")
plotDir      <- Sys.getenv("plotDir")
pdfDir       <- Sys.getenv("pdfDir")
centerWidth  <- as.numeric(Sys.getenv("centerWidth"))
minTLen      <- as.numeric(Sys.getenv("minTLen"))
maxTLen      <- as.numeric(Sys.getenv("maxTLen"))
minQual      <- as.numeric(Sys.getenv("minQual"))

# set variables
widthPadding <- (centerWidth - 1) / 2;
idx_col <- 'darkred'
oth_col <- 'blue'
non_col <- 'grey'
tlim    <- c(minTLen, maxTLen)
mlim    <- c(-15,15)
flim    <- c(0,1)
#dlim    <- tlim
dlim    <- c(0, 1)
vlim    <- tlim * 0.75

# get the data
write("  loading data", file=stderr())
d <- read.table(plotFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(d) <- c('chrom', 'region', 'group', 'grpId', 'cnvType', 'cnvSize', 'cnvFrac', 'idx_sample', 'idx_bin',
                 'end', 'sample', 'n_log_p', 'deltaCDF_sgn', 'nFrags');
s <- read.table(spanFile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
colnames(s) <- c('grpId', 'idx_sample', 'sample', 'start', 'end', 'tLen', 'qual');

# functions for CDFs (see mle.R for more explanation)
getDist <- function(sample){
    distFile        <- paste(pdfDir, "/svtools.pdf.", sample, ".all.gz", sep="")
    tmp             <- read.table(distFile, header=FALSE, sep="\t", stringsAsFactors=FALSE)
    colnames(tmp)   <- c('tLen', 'freq')
    tmp  
}
weightFreqs <- function(freq, tLen){
    wF <- freq * abs(tLen)
    wF / sum(wF)
}

# plot functions
addVLines <- function(){
    abline(v=c(left, right, center)/1e6, col="grey", lty=2)
}
plotDeltaCDF <- function(){
    legend("topleft", idx_sample, col=idx_col, lwd=1, bty="n", cex=1)
    abline(h=0, col="black", lty=1)    
    addVLines()
    for (sample in unique(ddd[,'sample'])){    
        dddd <- ddd[ddd$sample==sample,]
        if(nrow(dddd)>0){
            dddd[,'deltaCDF_sgn'] <- ifelse(dddd$deltaCDF_sgn < mlim[1], mlim[1], dddd$deltaCDF_sgn)        
            dddd[,'deltaCDF_sgn'] <- ifelse(dddd$deltaCDF_sgn > mlim[2], mlim[2], dddd$deltaCDF_sgn)
            if(sample==idx_sample) {
                #lines (dddd$end/1e6, dddd$deltaCDF_alt, col="black")
                #points(dddd$end/1e6, dddd$deltaCDF_alt, col="black", pch=20)    
                lines (dddd$end/1e6, dddd$deltaCDF_sgn, col=idx_col)
                points(dddd$end/1e6, dddd$deltaCDF_sgn, col=idx_col, pch=20)
            } else if(sample==oth_sample) {
                lines (dddd$end/1e6, dddd$deltaCDF_sgn, col=oth_col)   
            } else {
                lines (dddd$end/1e6, dddd$deltaCDF_sgn, col=non_col)
            }            
        }
    }  
}
plotSpans <- function(sample, col){
    legend("bottomleft", sample, col=col, lwd=1, bty="n", cex=1)
    abline(h=0, col="black", lty=1)    
    addVLines()    
    sss <- ss[ss$sample==sample,]
    if(nrow(sss)>0){
        col <- ifelse(sss$qual>=minQual, col, "darkgrey")
        segments(sss$start/1e6, sss$tLen, sss$end/1e6, sss$tLen, col=col, lwd=1, lend=2)           
    }
}
getSampleCDF <- function(sample, tL){
    dist      <- getDist(sample)    
    dist$freq <- weightFreqs(dist$freq, dist$tLen)
    sizeBin   <- dist[2,'tLen'] - dist[1,'tLen']
    tLens     <- seq(minTLen, maxTLen, sizeBin)
    pdf_ref   <- stepfun(dist$tLen, c(0, dist$freq))
    cum_sum   <- cumsum(dist$freq)
    cdf_ref   <- stepfun(dist$tLen, c(0, cum_sum))
    ref_inv   <- stepfun(cum_sum, c(dist[1,'tLen'], dist$tLen))    
    cdf_act   <- ecdf(tL)
    act_inv   <- stepfun(cdf_act(tLens), c(tLens[1], tLens))
    cnvType   <- ddd[1,'cnvType']
    cnvSize   <- ddd[1,'cnvSize']
    cnvFrac   <- ddd[1,'cnvFrac']
    pdf_alt   <- pdf_ref(tLens) * (1 - cnvFrac) + pdf_ref(tLens + cnvSize) * cnvFrac
    pdf_alt   <- pdf_alt / sum(pdf_alt)
    cdf_alt   <- stepfun(tLens, c(0, cumsum(pdf_alt)))
    list(tLens=tLens, cnvType=cnvType, cnvSize=cnvSize, cnvFrac=cnvFrac,
         cdf_ref=cdf_ref, ref_inv=ref_inv, cdf_act=cdf_act, act_inv=act_inv, cdf_alt=cdf_alt)
}
plotCDF_diffs <- function(){
    abline(h=0, col="black", lty=1) 
    abline(v=0, col="black", lty=1)  
    for (sample in unique(tttt[,'sample'])){
        ttttt <- tttt[tttt$sample==sample,]
        if(nrow(ttttt)>0){
            cdf   <- getSampleCDF(sample, ttttt$tLen)
            col   <- if(sample==idx_sample){
                idx_col
            } else if(sample==oth_sample) {
                oth_col
            } else {
                non_col
            }
            frqs <- seq(0, 1, 0.01)
            dltTLens <- cdf[['act_inv']](frqs) - cdf[['ref_inv']](frqs)       
            lines(dltTLens, frqs, col=col)
            #frq_act <- cdf[['cdf_act']](cdf[['tLens']])
            #tls <- ifelse(frq>=0&frq<=1, cdf[['ref_inv']](frq), NA)
            #lines(cdf[['tLens']], -(cdf[['tLens']]-tls), col=col)            
        }
    }
}
plotCDFs <- function(sample, col){
    abline(h=0, col="black", lty=1) 
    abline(v=0, col="black", lty=1)    
    ttttt <- tttt[tttt$sample==sample,]
    if(nrow(ttttt)>0){
        cdf   <- getSampleCDF(sample, ttttt$tLen)
        tL  <- cdf[['tLens']]
        ref <- cdf[['cdf_ref']](tL)
        act <- cdf[['cdf_act']](tL)
        polygon(c(tL,rev(tL)),c(ref,rev(act)),col="lightgrey",lty=0)    
        if(sample==idx_sample){
            if(cdf[['cnvType']]=='del'){cdf[['cnvFrac']] <- 1-cdf[['cnvFrac']]}
            abline(h=cdf[['cnvFrac']], col="grey", lty=2) 
            abline(v=-as.numeric(cdf[['cnvSize']]), col="grey", lty=2)
            alt <- cdf[['cdf_alt']](tL)
            #polygon(c(tL,rev(tL)),c(alt,rev(act)),col="lightgreen",lty=0)
        } 
        lines(tL, ref, col="black")
        lines(tL, act, col=col)    
        if(sample==idx_sample){
            lines(tL, alt, col="darkgreen")  
        }         
    }
}

# make the plots
for (grpId in unique(d$grpId)){
    
    # get this index event
    write(paste("  ", grpId, sep=""), file=stderr())    
    dd      <- d[d$grpId==grpId,]
    ss      <- s[s$grpId==grpId,]
    group   <- dd[1,'group']    
    
    # set the coordinates
    chrom   <- dd[1,'chrom']
    start   <- min(dd[,'end'])
    end     <- max(dd[,'end'])
    padding <- (end - start) / centerWidth    
    left    <- start + padding * widthPadding
    right   <- end   - padding * widthPadding
    center  <- left + (right - left) / 2
    clim    <- c(start, end) /  1e6
    idx_bin <- dd[1,'idx_bin']
    ttt     <- ss[ss$start<idx_bin&ss$end>idx_bin&ss$qual>=minQual,]
    
    # determine the index and other samples
    idx_samples <- unique(dd$idx_sample) # can be more than one if region has multiple candidate samples
    idx_sample  <- idx_samples[1]
    ddd <- dd[dd$idx_sample==idx_sample,] # only need to plot the trace from one call
    tttt <- ttt[ttt$idx_sample==idx_sample,]
    oth <- ddd[ddd$sample!=idx_sample&ddd$end>=left&ddd$end<=right,]
    dlt <- abs(oth[,'deltaCDF_sgn'])
    oth_sample <- oth[which(dlt==max(dlt))[1],'sample']

    # open the jpg file
    jpgFile <- paste(ddd[1,'region'], group, grpId, "jpg", sep=".")
    jpgFile <- paste(plotDir, jpgFile, sep="/")        
    jpeg(jpgFile,
         width=5, height=5, units="in", pointsize=9,
         res=600, quality=100)    
    
    # set the composite layout
    layout(matrix(c(1,4,2,5,3,6), 3, 2, byrow = TRUE), width=c(0.7,0.3))

    # plot the MLE p-value traces
    par(mar=c(0.1, 4.1, 4.1, 0.5))
    plot(0, 0, type="n", xlim=clim, ylim=mlim,
         xaxt="n",
         ylab="deltaCDF",
         main=paste(paste(group, " #", grpId, ", ",
                          ddd[1,'cnvType'], " ", as.numeric(ddd[1,'cnvSize']), " bp", sep=""))
    )
    plotDeltaCDF()

    # plot the index sample spans
    par(mar=c(2.05, 4.1, 2.05, 0.5))
    plot(0, 0, type="n", xlim=clim, ylim=tlim,
         xaxt="n",
         ylab="TLEN"
    )
    plotSpans(idx_sample, idx_col)
    
    # plot the other/control sample spans
    par(mar=c(4.1, 4.1, 0.1, 0.5))
    plot(0, 0, type="n", xlim=clim, ylim=tlim,
         xlab=paste(chrom, "coordinate (Mb)", sep=" "),
         ylab="TLEN"
    )
    plotSpans(oth_sample, oth_col)
    
    # plot the CDF differences for all samples
    par(mar=c(0.1, 4.1, 4.1, 0.5))
    plot(0, 0, type="n", xlim=vlim, ylim=dlim,
         xaxt="n",
         ylab="cumulative fraction"
    ) 
    plotCDF_diffs()
      
    # plot the index sample CDFs
    par(mar=c(2.05, 4.1, 2.05, 0.5))
    plot(0, 0, type="n", xlim=vlim, ylim=flim,
         xaxt="n",
         ylab="cumulative fraction"
    ) 
    plotCDFs(idx_sample, idx_col)
    
    # plot the other/control sample CDFs
    par(mar=c(4.1, 4.1, 0.1, 0.5))
    plot(0, 0, type="n", xlim=vlim, ylim=flim,
         xlab="TLEN or delta(TLEN)",
         ylab="cumulative fraction"
    )
    plotCDFs(oth_sample, oth_col)

    # finish the jpg
    graphics.off()
}








#plotPValues <- function(){
#    legend("topleft", idx_sample, col=idx_col, lwd=1, bty="n", cex=1)
#    abline(h=0, col="black", lty=1)    
#    abline(h=c(3,5), col="grey", lty=2)
#    addVLines()
#    for (sample in unique(ddd[,'sample'])){    
#        dddd <- ddd[ddd$sample==sample,]
#        if(sample==idx_sample) {
#            lines (dddd$end/1e6, dddd$nlog_ks_p_alt, col="black")
#            points(dddd$end/1e6, dddd$nlog_ks_p_alt, col="black", pch=20)            
#            lines (dddd$end/1e6, dddd$nlog_ks_p_ref, col=idx_col)
#            points(dddd$end/1e6, dddd$nlog_ks_p_ref, col=idx_col, pch=20)
#        } else if(sample==oth_sample) {
#            lines (dddd$end/1e6, dddd$nlog_ks_p_ref, col=oth_col)   
#        } else {
#            lines (dddd$end/1e6, dddd$nlog_ks_p_ref, col=non_col)
#        }
#    }  
#}