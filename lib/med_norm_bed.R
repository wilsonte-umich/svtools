
d <- read.table("stdin", header=FALSE, sep="\t")
for(chr in unique(d$V1)){
    dd <- d[d$V1==chr,]
    m <- median(dd$V5)
    dd$V7 <- dd$V5 / m
    write.table(dd, file="", quote=FALSE, sep="\t", row.names = FALSE, col.names = FALSE, append=TRUE)    
}

## find coverage mode
#ux <- unique(d5)
#m <- ux[which.max(tabulate(match(d5, ux)))]
#cat(m, " = mode\n")
#
#cat(median(d5), " = median\n")
#
## generate a coverage histogram
#h <- hist(d5, max(d5), plot=FALSE)
#freqs <- h$counts/nrow(d)
#
## curve fit a gaussian
#a <- max(freqs)
#m <- sum(freqs * h$mids)
#freqSum <- cumsum(freqs)
#s <- min(h$mids[freqSum>=0.841]) - m
#guassian <- freqs ~ a*(exp(-((mids-m)^2)/(2*(s^2))))
#curveFit <- try(nls(guassian,
#                    data=data.frame(freqs=freqs,mids=h$mids),
#                    start=list(a=a,m=m,s=s)),
#                silent=TRUE)
#curveFitOK <- !inherits(curveFit, "try-error")
#if(curveFitOK){
#    coef <- coef(curveFit)
#    a <- coef[['a']]
#    m <- coef[['m']]
#    s <- coef[['s']]
#} else {
#    cat("ERROR: curve fit failed", "\n", file=stderr())
#}
#
#save(h,file="test.data")
#cat(a, m, s, "\n")
