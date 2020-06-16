#
## get passed arguments
#fracAlt  <- as.numeric(Sys.getenv("fracAlt"))
#
## does this need lambda?
#
## set variables
#zeroProb <- 0.001 # effectively the base error rate in sequencing
#
## set columns
#nRef_col    <- 7
#nAlt_col    <- 8
#type_col    <- 9
#fracAlt_col <- 10
#
## calculate probablity of simple het vs. alternative model
#probs <- function(v){
#    nRef   <- as.numeric(v[nRef_col])
#    nAlt   <- as.numeric(v[nAlt_col])    
#    nReads <- nRef + nAlt     
#    
#    res <- if(nReads > 0) {     
#        type        <- v[type_col]
#        fracAlt     <- as.numeric(v[fracAlt_col])
#        fracAlt_het <- 1/2 * (fracAlt / 0.5)
#        pHet        <- -log10(dbinom(nAlt, nReads, fracAlt_het))        
#        
#        
#        
#
#        pCnv <- if(type == 'refDup'){
#            fracAlt_refDup <- 1/3 * (fracAlt / 0.5)
#            -log10(dbinom(nAlt, nReads, fracAlt_refDup))
#        } else if(type == 'altDup'){
#            fracAlt_altDup <- 2/3 * (fracAlt / 0.5)
#            -log10(dbinom(nAlt, nReads, fracAlt_altDup))
#            
#            
#        } else if(type == 'refDel'){
#            -log10(zeroProb ** nRef)
#        } else if(type == 'altDel'){
#            -log10(zeroProb ** nAlt)
#        }
#        
#        c(pHet, pCnv, pHet-pCnv)
#    } else {
#        c(0,0,0)
#    }
#    paste(res, collapse="\t")
#    
#    # deal with Inf
#    
##chr19   3111888 4436682 1       741     .       613     476     het     2.10082869594257        Inf     -Inf
##chr19   4436755 4438709 2       9       .       0       10      refDel  3.37527763618577        0       3.37527763618577
##chr19   4438834 4642166 3       232     .       209     167     het     1.46284092248514        Inf     -Inf
##chr19   4642181 4642613 4       12      .       10      0       altDel  2.67363103339128        0       2.67363103339128
#}
#
## get data
#data <- read.table(file="stdin", header=FALSE, sep="\t")
#data[,fracAlt_col] <- apply(data, 1, probs)
#
## calculate snp-specific emission probabilities
#write.table(data, file="", quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
