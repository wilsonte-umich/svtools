
# get passed arguments
args <- commandArgs(TRUE)
tmpFile <- args[1]
covFile <- args[2]
medTLen <- as.numeric(args[3])
chromMax <- as.numeric(args[4])
chrom <- args[5]

# open input stream
zcat <- paste('zcat', tmpFile)
awk <- paste("awk '", '$1=="', chrom, '"', "'", sep="")
cut <- 'cut -f 2,3'
inPipe <- paste(zcat, awk, cut, sep=" | ")
con <- pipe(inPipe, open="r")

# create the coverage map
nBases  <- rep(0, chromMax)
while(length(v <- scan(con,what=integer(),nlines=1,quiet=TRUE))>0){
    is <- v[1]:v[2]
    nBases[is] <- nBases[is] + 1;
}

# close stream
close(con)

# normalize counts to fragment length and save R vector to file
nBases <- round(nBases[1:chromMax] / medTLen, 4)
save(nBases, file=covFile, compress=TRUE)

cat("OK\n")
