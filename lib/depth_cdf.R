
# get passed arguments
tmpFile <- Sys.getenv("tmpFile")
rFile   <- Sys.getenv("rFile")
jpgFile <- Sys.getenv("jpgFile")

# get the data
d <- read.table(tmpFile, header=TRUE, sep="\t")

# determine data limits
median <- median(d[,1])
xlim   <- c(median - median/2, median + median/2)

# make the plot
jpeg(jpgFile,
     width=4, height=4, units="in", pointsize=10,
     res=600, quality=100)
plot(0, 0, type="n",
     xlim=xlim, ylim=c(0,1),
     xlab="Bin Count",
     ylab="Cumulative Fraction")
abline(v=median, col="grey", lty=2)
abline(v=median-median/5, col="grey", lty=2)
abline(v=median+median/5, col="grey", lty=2)
abline(h=0.2, col="grey", lty=2)
abline(h=0.5, col="grey", lty=2)
abline(h=0.8, col="grey", lty=2)
col <- 1
cols <- numeric()
for (i in 1:ncol(d)){
    write(paste(colnames(d)[i], round(sd(d[,i]),3), sep=" => sd "), file=stderr())
    counts <- round(d[,i], 0)
    xy     <- aggregate(d[,i], list(counts), length)
    xy[,2] <- cumsum(xy[[2]] / sum(xy[[2]]))
    lines(xy[[1]], xy[[2]], col=col)
    cols <- c(cols, col)
    col <- col + 1
}
legend("bottomright", colnames(d),
       col=cols, lwd=1, bty="n", cex=0.8)

