
# get passed arguments
dataFile  <- Sys.getenv("dataFile")
#lambda    <- as.numeric(Sys.getenv("lambda"))
#fracAlt   <- as.numeric(Sys.getenv("fracAlt"))
errorRate <- as.numeric(Sys.getenv("errorRate"))

# set variables
nReadCorr      <- c(2/2, 3/2, 3/2, 1/2, 1/2)
pFracAlt_zero  <- rep(1, 5)
#fracAlt_het    <- 1/2 * (fracAlt / 0.5)
#fracAlt_refDup <- 1/3 * (fracAlt / 0.5)
#fracAlt_altDup <- 2/3 * (fracAlt / 0.5)

# calculate eProbs as states het, refDup, altDup, refDel, altDel
eProbs <- function(v){
    nReads   <- v[1] + v[2]
    pNReads  <- dpois(nReads, v[3] * nReadCorr)
    pFracAlt <- if(nReads == 0) {
        pFracAlt_zero        
    } else {    
        c(dbinom(v[2], nReads, 1/2 * (v[4] / 0.5)),
          dbinom(v[2], nReads, 1/3 * (v[4] / 0.5)),
          dbinom(v[2], nReads, 2/3 * (v[4] / 0.5)),
          errorRate ** v[1],
          errorRate ** v[2])         
    }
    ep <- pNReads * pFracAlt
    paste(c(paste(ep, collapse=","), paste(log10(ep), collapse=",")), collapse="\t")
}

# get data (eProb, i, nRef, nAlt, lambda, fracAlt)
data <- read.table(dataFile, header=FALSE, sep="\t")
data[,1] <- apply(data[,c(3:6)], 1, eProbs)

# calculate snp-specific emission probabilities
write.table(data, file="", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
