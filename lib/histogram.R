
# get passed arguments
args <- commandArgs(TRUE)
imgFile <- args[1]
sample <- args[2]
maxX <- as.numeric(args[3])
xlab <- args[4]

# load data
d <- read.table("stdin", header=FALSE, sep="\t");
n <- sum(d[,2])
freq <- d[,2] / n
maxF <- max(freq, na.rm=TRUE)
jpeg(imgFile,
     width=2.5, height=2.5, units="in", pointsize=9,
     res=600, quality=100)
plot(d[,1], freq,
     type="h", main=sample,
     xlim=c(0,maxX), xlab=xlab,
     ylim=c(0,maxF*1.05), ylab="Frequency")
