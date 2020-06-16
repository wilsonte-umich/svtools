
# passed arguments
args         <- commandArgs(trailingOnly=TRUE)
nCnvLeadCol  <- as.numeric(args[1])

# load data
d        <- read.table(file('stdin'), header=FALSE, sep="\t", stringsAsFactors=FALSE)
maxCol   <- ncol(d)
nSamples <- (maxCol - nCnvLeadCol) / 2
minCol   <- nCnvLeadCol + nSamples + 1

# calculate probability ratio of gain/loss vs. neutral for each sample
getPRatio <- function(v){ 
    cnvDcn  <- as.numeric(strsplit(v[2], ',')[[1]])       
    if(length(cnvDcn) >= 2){
        altDcn <- if(as.numeric(v[1]) > 0) { 1 } else { -1 }
        tryCatch({
            round(log10(t.test(cnvDcn, mu=altDcn)$p.value / t.test(cnvDcn, mu=0)$p.value), 3)
        }, error = function(e) {
            0
        })
    } else {
        0
    }
}
for (col in minCol:maxCol){
    d[,col] <- apply(d[,c(col-nSamples,col)], 1, getPRatio)
}

# return data on stdout
write.table(d, file=stdout(), quote=FALSE, sep = "\t", row.names=FALSE, col.names=FALSE)
