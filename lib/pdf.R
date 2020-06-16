
# get passed arguments
sample      <- Sys.getenv("sample")
libraries   <- strsplit(Sys.getenv("libraries"), ",")[[1]]
distFiles   <- strsplit(Sys.getenv("distFiles"), ",")[[1]]
jpgFile     <- Sys.getenv("jpgFile")

# get the data
libI <- 1:length(libraries)
dist <- list()
maxF <- 0
for (i in libI){
    dist[[i]] <- read.table(distFiles[i], header=FALSE, sep="\t", stringsAsFactors=FALSE)
    maxF <- max(maxF, dist[[i]][,2])
}

# get 50th, 95th, and 99th percentile lines
cum    <- cumsum(dist[[1]][,2])
tLen50 <- dist[[1]][max(which(cum<=0.5)),1]
tLen95 <- dist[[1]][max(which(cum<=0.95)),1]
tLen99 <- dist[[1]][max(which(cum<=0.99)),1]

# make the plot
jpeg(jpgFile,
     width=4, height=4, units="in", pointsize=10,
     res=600, quality=100)
plot(dist[[1]][,1], dist[[1]][,2], type="n",
     ylim=c(0,maxF),
     xlab="BAM TLEN",
     ylab="Fraction of FR Pairs",
     main=sample)
abline(v=tLen50, lty=2)
abline(v=tLen95, lty=2)
abline(v=tLen99, lty=2)
col  <- 1
cols <- numeric()
for (i in libI){
    lines(dist[[i]][,1], dist[[i]][,2], col=col)
    cols <- c(cols, col)
    col <- col + 1
}
legend("topright", libraries, col=cols, lwd=1, bty="n", cex=0.8)




