
# reporting functions
commify <- function(x) { # pretty comma-delimited formatting of numbers
    format(x, big.mark=",", scientific=FALSE)
}
commifyInt <- function(x) {
    commify(round(x, 0))
}
reportCount <- function(count, msg){
    message(paste(commifyInt(count), msg, sep=" "))
}
head.table <- function(table, nrow=10, ncol=10){
    print(table[1:nrow,1:ncol])
}

# file functions
utility <- "svtools"
get.file <- function(dir, type, name, suffix, fileCommand=NULL){
    if(is.null(fileCommand)) fileCommand <- command # default to current action
    prefix <- paste(dir, "/", utility, ".", fileCommand, sep="")
    paste(prefix, type, name, suffix, sep=".") 
}
get.inFile <- function(type, name, suffix, fileCommand=NULL){
    get.file(args$IN_DIR, type, name, suffix, fileCommand)
}
get.outFile <- function(type, name, suffix, fileCommand=NULL){
    get.file(args$OUT_DIR, type, name, suffix, fileCommand)
}
get.plotFile <- function(type, name){
    get.file(args$PLOT_DIR, type, name, "png")
}
read.bgz <- function(type, name, fileCommand=NULL, header=FALSE){
    file <- get.inFile(type, name, "bgz", fileCommand)
    message(paste("reading file:", file, sep=" ")) 
    d <<- read.table(file, header=header, sep="\t",
                     stringsAsFactors=FALSE, comment.char="",)
    head.table(d)
}
write.bgz <- function(table, type, name, col.names=FALSE){
    file <- get.outFile(type, name, "bgz")
    message(paste("writing file:", file, sep=" "))    
    bgz <- paste("bgzip -c >", file, sep=" ")    
    write.table(table, file=pipe(bgz), quote=FALSE, sep="\t",
                row.names=FALSE, col.names=col.names)
    tabix <- paste("tabix -p bed", file, sep=" ")
    system(tabix)
    head.table(table)
}

# vector functions
collapseVector <- function(v, n) { # sum every n adjacent elements of a vector 
    unname(tapply(v, (seq_along(v)-1) %/% n, sum))
}

# debugging actions
DEBUG_FILE <- get.outFile("TEST", "TEST", "RData")
debug.file <- function(action="r"){
    if(action == "r"){
        load(DEBUG_FILE, envir=topenv())
        print(objects(envir=topenv()))
    } else if(action == "w") {
        save.image(DEBUG_FILE)
        quit("no")
    }
}

# plotting functions
plotDim   <- 900
pointsize <- 12
plotHistogram <- function(v, type, xlab,
                          xlim=range(v), vLine=NULL, breaks=100){
    p <- v[v >= xlim[1] & v <= xlim[2]]
    file <- get.plotFile(type, args$SAMPLE)
    png(file, width=plotDim, height=plotDim,
        units="px", pointsize=pointsize, type = c("cairo"))    
    hist(p, breaks=breaks, main="", xlab=xlab, xlim=xlim, ylab="Frequency")
    if(!is.null(vLine)) abline(v=vLine, col="red")
    graphics.off()
}

# chromosome functions
getOrderedChroms <- function(chroms, includeY=TRUE){
  chroms_ <- sub("chr", "",  unique(chroms))
  isAutosome <- !is.na(suppressWarnings(as.numeric(chroms_)))
  chroms <- paste("chr", c(sort(as.numeric(chroms_[isAutosome])), sort(chroms_[!isAutosome])), sep="")
  if(!includeY) chroms <- chroms[!(chroms %in% "chrY")]
  chroms
}

#> m
#     [,1] [,2] [,3]
#[1,]    1    3    5
#[2,]    2    4    6
#> t(apply(m, 1, '*', 2:4)) <<<<< multiply rows by a vector
#     [,1] [,2] [,3]
#[1,]    2    9   20
#[2,]    4   12   24
#> apply(m, 2, '*', 2:3) <<<<< multiply columns by a vector
#     [,1] [,2] [,3]
#[1,]    2    6   10
#[2,]    6   12   18


