
# this script find runs of usable probes that predict structural variants
# in one sample relative to an input model built on many samples

# runs of probes may be CNVs, LOH, or both

# function for final printing of SVs BED file
printSVs <- function(){
    message("writing SVs BED")
    header <- paste(SV_COL, collapse="\t")
    header <- paste("#", header, sep="") # comment the header for bgzip
    write(header, file=SVS_FILE)         # print it
    pipe <- paste("awk 'BEGIN{OFS=\"\\t\"}{$2+=0;$3+=0;print $0}' >> ", SVS_FILE, sep="")    
    nSVOut <- nrow(SVs)  
    if(nSVOut>0){     
        SVs$SAMPLE <- MODEL_NAME # rename a couple of columns for header
        SVs$SV_ID  <- SVs$segN
        for(col in c('LRR_MEAN','LRR_SD','ZYG_MEAN','ZYG_SD')){
            SVs[,col] <- round(SVs[,col],3) # pretty-round some numbers
        }    
        write.table(                         # then append the data rows
            SVs[,SV_COL],
            file=pipe(pipe),
            quote=FALSE,
            sep="\t",
            row.names=FALSE,
            col.names=FALSE
        )    
    }
    message(paste("   ", nSVOut, "structural variants written to disk"))
    quit("no")
}

# set parameters
MIN_PROBES_SV  <- 5   # require this many probes in all sub-runs (typically adjacent)
MIN_INF_BITS   <- 7   # require this many bits of information, where a bit = 1 LRR OR l informative ZYG
MIN_REL_LL     <- 10  # reject probes that don't show at least this relative log-likelihood
MAX_GAP_BP     <- 1e5 # merge runs with less than this many intervening non-anomalous bp
N_FLANK_PROBES <- 10  # number of probes on either side of SV span used in future to determine endpoint uniqueness

# prepare row information for finding SV runs
message("identifying anomalous probes")
prbI <- 0 # a sequential numerical probe identifier for usable probes
di   <- sapply(d$IS_NA==1, function(IS_NA){ if(IS_NA) NA else prbI <<- prbI + 1 })
by   <- list(d$CHROM)
maxDI  <- max(di, na.rm=TRUE)                # get limits on probe identifiers
minDIC <- aggregate(di, by, min, na.rm=TRUE) # here by chrom
maxDIC <- aggregate(di, by, max, na.rm=TRUE)
colnames(minDIC) <- c('CHROM','i')
colnames(maxDIC) <- c('CHROM','i')
ducp <- d[d$IS_NA==0,c('CHROM','POS')]       # coordinates of usable probes
cnc  <- d$CNC != 0                           # probes with a copy number change
loh  <- astasms[d$AS_OUT] < astasms[d$AS_IN] # probes with LOH, even if copy number neutral
anom <- d$IS_NA==0 & (cnc | loh)  # usable probes with any called deviation from the input model
rm(cnc, loh)                      # TODO: could only report CNC and ignore LOH

# identify finest granularity runs of contiguous anomalous probes
message("finding runs of anomalous probes")
dd     <- d [anom,]
ddi    <- di[anom] # row numbers from d for rows present in dd (NOT row numbers from dd)
segN_  <- 0 # output run identifier
prv    <- 1:(nrow(dd)-1) # set up probe adjacency comparisons
nxt    <- 2: nrow(dd)    
diff_chrom <- c(TRUE, dd[nxt,'CHROM'] != dd[prv,'CHROM'])
any_gap    <- c(TRUE, ddi[nxt] - ddi[prv] > 1)
new_cn     <- c(TRUE, dd[nxt,'CN_OUT'] != dd[prv,'CN_OUT']) # will therefore fuse CN2_11 and CN2_02, e.g.
dd$segN <- sapply( diff_chrom | any_gap | new_cn,
    function(inc){ # break on new chromosome, or _any_ probe gap, or if CN changes
        if(inc) segN_ <<- segN_ + 1
        segN_
    })
rm(prv, nxt, diff_chrom, any_gap, new_cn)

# filter away trivially short probe runs
# also, filter away SVs claimed as CNN LOH where no probes are informative
# such regions result from the vagaries of HMM assignment within a CN run
message(paste("limiting output to runs with anomalous probe count >=", MIN_PROBES_SV))
SVs <- data.frame(segN=sort(unique(dd$segN)))
by  <- list(segN=dd$segN)
SVs$N_PROBES <- aggregate(1:nrow(dd), by, function(x){
    length(which( 
        dd[x,'CNC'] != 0 | (dd[x,'INF_IN']  == 'yes' &
                            dd[x,'INF_OUT'] == 'no')
    ))
})[[2]]
SVs$N_INF <- aggregate(1:nrow(dd), by, function(x){
    length(which(dd[x,'INF_IN']  == 'yes'))
})[[2]]
keep_SVs <- function(flt){ # function for filtering SVs
    SVs <<- SVs[flt,]
    if(nrow(SVs)==0){ # abort if no SVs remain to avoid subsequent errors
        message("    no structural variants remain after filtering")
        printSVs()
    } else { # otherwise, continue filtering and proceed
        flt  <- dd$segN %in% SVs$segN
        dd  <<- dd [flt,]
        ddi <<- ddi[flt]
        by  <<- list(segN=dd$segN)            
    }
}
keep_SVs(SVs$N_PROBES >= MIN_PROBES_SV &           # give more weight to informative probes
         SVs$N_PROBES + SVs$N_INF >= MIN_INF_BITS) # i.e. allow shorter runs if contains informative

# calculate likelihoods for each SV run
# measure of SV strength, degree of anomaly rel. to array
message("calculating SV likelihoods")
SVs[,c('MDL_LL','SV_LL')] <- t(sapply(1:nrow(SVs), function(i){
    ddi <- which(dd$segN %in% SVs[i,'segN']) # row numbers from dd
    obs <- dd[ddi,c('ot','os')] # observation indices  
    sapply(c(
        HMM_likelihood_obs(HMM, obs, dd[ddi,'AS_IN'] ), # likelihood of input model path
        HMM_likelihood_obs(HMM, obs, dd[ddi,'hs'])      # likelihood of SV path
    ), round, 1)
}))
SVs$REL_LL <- SVs$SV_LL - SVs$MDL_LL  

# filter away probe runs without sufficiently greater likelihood relative to input model
message(paste("limiting output to runs with relative log likelihood >=", MIN_REL_LL))
keep_SVs(SVs$REL_LL >= MIN_REL_LL)

# establish the coordinate spans of the SV runs
message("aggregating SV parameters")
FIRST <- function(x) x[1]
SVs$CHROM  <- aggregate(dd$CHROM,  by, FIRST)[[2]]
SVs$START  <- aggregate(dd$START,  by, min)[[2]]
SVs$END    <- aggregate(dd$POS,    by, max)[[2]]
SVs$SPAN   <- SVs$END - SVs$START
SVs$STRAND <- "+"

# set the SV flanking endpoints in bp corresponding to N_FLANK_PROBES
# used later for uniqueness determination
minFlnkI <- sapply(aggregate(ddi-N_FLANK_PROBES, by, min)[[2]], max, 1)
maxFlnkI <- sapply(aggregate(ddi+N_FLANK_PROBES, by, max)[[2]], min, maxDI)
minFlnkI <- sapply(1:nrow(SVs), function(j){ # be sure not to overrun the chromsome
    i <- minFlnkI[j]
    c <- SVs[j,'CHROM']
    if(ducp[i,'CHROM'] != c){ minDIC[minDIC$CHROM==c,'i'] } else { i } 
})
maxFlnkI <- sapply(1:nrow(SVs), function(j){
    i <- maxFlnkI[j]
    c <- SVs[j,'CHROM']
    if(ducp[i,'CHROM'] != c){ maxDIC[maxDIC$CHROM==c,'i'] } else { i } 
})
SVs$FLANK_START <- ducp[minFlnkI,'POS'] 
SVs$FLANK_END   <- ducp[maxFlnkI,'POS']

# determine ~fixed SV properties of each run  
COLLAPSE <- function(x) paste(sort(unique(x)), collapse=",")
CONTAINS <- function(x) ifelse('yes' %in% x, 'yes', 'no')
ANY      <- function(x) ifelse(any(x), 'yes', 'no')
SVs$AS_IN  <- aggregate(dd$AS_IN,  by, COLLAPSE)[[2]]
SVs$AS_OUT <- aggregate(dd$AS_OUT, by, COLLAPSE)[[2]]
SVs$CN_IN  <- aggregate(dd$CN_IN,  by, COLLAPSE)[[2]]
SVs$CN_OUT <- aggregate(dd$CN_OUT, by, COLLAPSE)[[2]]
SVs$CNC    <- aggregate(dd$CNC,    by, COLLAPSE)[[2]]
SVs$LOH    <- aggregate(astasms[dd$AS_OUT] < astasms[dd$AS_IN], by, ANY)[[2]]

# aggregate probe properties over the SV run
# only aggregate zygosity for probes that were informative in the input model
ZYG <- ifelse(dd$INF_IN=='yes', dd$ZYG, NA)
SVs$LRR_MEAN <- aggregate(dd$LRR,    by, mean)[[2]]
SVs$LRR_SD   <- aggregate(dd$LRR,    by, sd)[[2]]
SVs$ZYG_MEAN <- aggregate(ZYG, by, mean, na.rm=TRUE)[[2]]
SVs$ZYG_SD   <- aggregate(ZYG, by, sd,   na.rm=TRUE)[[2]]
SVs[is.nan(SVs$ZYG_MEAN),'ZYG_MEAN'] <- NA
SVs[is.nan(SVs$ZYG_SD),  'ZYG_SD']   <- NA
rm(ZYG)

# calculate inferred (not called) copy numbers for each SV
calc_cn <- function(LRR) 2 * 2**LRR # expects adjusted log2 LRR values
CN_CLC <- aggregate(dd$LRR, by, function(LRR){
    paste(calc_cn(LRR), collapse=",")
})[[2]]
CN_corr <- list()
for (CN in unique(d$CN_IN[!is.na(d$CN_IN)])){
    CN_ <- as.character(CN)
    CN_corr[[CN_]] <- CN / median(calc_cn(d[d$IS_NA==0 &
                                            d$CN_IN==d$CN_OUT &
                                            d$CN_IN==CN,'LRR']), na.rm=TRUE)
}
SVs$CN_CLC <- 0.0
SVs$CN_P   <- 0.0
for(i in 1:nrow(SVs)){
    cn_out  <- as.numeric(SVs[i,'CN_OUT'])
    cn_out_ <- as.character(cn_out)      
    cn_corr <- if(is.null(CN_corr[[cn_out_]])) 1 else CN_corr[[cn_out_]]  
    cn_clc <- as.numeric(strsplit(CN_CLC[i], ",")[[1]]) * cn_corr
    SVs[i,'CN_CLC'] <- round(median(cn_clc, na.rm=TRUE), 2)
    #SVs[i,'CN_P']   <- round(wilcox.test(cn_clc, mu=cn_out)$p.value, 4)
}

# commit the final output file, which will have some SVs
printSVs()



## collect parameters
#LIB_DIR       <- '/home/wilsonte_lab/club_house/etc/garage/bin/wilson_ware/msvtools-0.0.1/lib'
##MODEL_NAME    <- 'U2OS1-12APH'
#MODEL_NAME    <- 'Scr1-1'
#SAMPLES_DIR   <- '/home/wilsonte_lab/club_house/data/Glover/Irene/samples'
#DATAFILE      <- paste(SAMPLES_DIR, '/', 'msvtools.segment.data.', MODEL_NAME, '.gz', sep="")
#PLOT_DIR      <- paste(SAMPLES_DIR, '/plots/', MODEL_NAME, sep="")
#PROBES_FILE   <- paste(SAMPLES_DIR, '/msvtools.segment.probes.', MODEL_NAME, '.bed.bgz', sep="")
#SEGMENTS_FILE <- paste(SAMPLES_DIR, '/msvtools.segment.segments.', MODEL_NAME, '.bed.bgz', sep="")
#SVS_FILE      <- paste(SAMPLES_DIR, '/msvtools.segment.SVs.', MODEL_NAME, '.bed', sep="")
#TRAIN_RDATA   <- '__NULL__'
#PLOT_PREFIX   <- paste(PLOT_DIR, MODEL_NAME, sep="/")
#
## save data for later development
#if(exists('saveRData')){
#    remove('saveRData')
#    RD_FILE <- paste(PROBES_FILE, "RData", sep=".")
#    save.image(file=RD_FILE)
#    quit("no") 
#    
## load data when required/requested
#} else if(!exists('d') | exists('forceLoad')){
#    RD_FILE <- paste(PROBES_FILE, '.RData', sep="")    
#    message(paste("loading data:", MODEL_NAME))
#    load(file=RD_FILE)
#    remove('saveRData')
#}


#logL <- function(ddi, AS_col){
#    sum( log( sapply(ddi, function(i) dd[i,which(colnames(dd)==dd[i,AS_col])] ) ) )
#}

## calculate probabilities for each SV run 
#message("calculating SV probabilities")
#lr.test <- function(i, nll_col, alt_col){ # https://cran.r-project.org/web/packages/extRemes
#    round(log(pchisq(2 * (SVs[i,alt_col] - SVs[i,nll_col]),
#                     df=SVs[i,'N_PROBES'],
#                     lower.tail=FALSE)), 3)  
#}
#SVs[,c('SV_P','MDL_P')] <- t(sapply(1:nrow(SVs), function(i){
#    c(lr.test(i, 'FWD_LL', 'SV_LL'),
#      lr.test(i, 'FWD_LL', 'MDL_LL'))
#}))
##null = the model with fewer parameters = FORWARD
##alt  = the model with more parameters, has the same or greater log-likelihood = SV PATH
##df   = n probes
