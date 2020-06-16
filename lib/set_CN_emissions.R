
# get passed arguments
refPloidy   <- as.numeric(Sys.getenv("refPloidy"))
maxCN       <- as.numeric(Sys.getenv("maxCN"))
eProbFile   <- Sys.getenv("eProbFile_CN")
refMean     <- as.numeric(Sys.getenv("refMean"))
refStdev    <- as.numeric(Sys.getenv("refStdev"))
binWidth    <- as.numeric(Sys.getenv("binWidth"))
maxBinSize  <- as.numeric(Sys.getenv("maxBinSize"))
minBinSize  <- as.numeric(Sys.getenv("minBinSize"))

# generate emission states = nCN + two limit states
CNs   <- c(1:maxCN) 
nCN   <- length(CNs)
nEmmStates <- nCN + 2

# generate observation states = nBins plus two limit states
bins  <- seq(minBinSize, maxBinSize, binWidth)
nBins <- length(bins)
nObsStates <- nBins + 2

# generate the emission probability matrix within HMM model boundaries
matrix <- matrix(0, nrow=nEmmStates, ncol=nObsStates)
for(CN in CNs){
    sc <- refPloidy / CN
    matrix[CN+1,] <- c(0, dnorm(bins, refMean * sc, refStdev * sc), 0) * binWidth 
}

# add the limit states
matrix[1,nObsStates] <- 1 
matrix[nEmmStates,1] <- 1

# print the results to file
write(t(matrix), file=eProbFile, ncolumns=nObsStates, append=FALSE, sep="\t")
