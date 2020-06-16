
# get passed agruments
createArgs <- as.list(Sys.getenv(c(
    "LIB_DIR",
    "HMM_FILE",
    "BIN_WEIGHT"
)))
createArgs$BIN_WEIGHT <- as.numeric(createArgs$BIN_WEIGHT)

# check for something to do
if(file.exists(createArgs$HMM_FILE)) {
    message(paste(createArgs$HMM_FILE, "already exists, nothing to do"))
    q("no")
}

# set state boundaries
message("setting model limits")
minCN   <- 1 # variable-width bins ALWAYS have a copy number (CN=0 spans should be subtracted as gaps)
modalCN <- 2 # assume base genome is diploid, e.g. human
maxCN   <- 4 # detect up to a doubling of baseline content in a region (anything higher reduced to 4.0)
CNs     <- minCN:maxCN
minN    <- 0
maxN    <- createArgs$BIN_WEIGHT * maxCN * 10 # a generous count range up to the highest allowed copy number
medians <- minN:(maxN*10)               # integer expected count states every 1/10th N for more HMM precision
Ns      <- minN:maxN                    # observed read count, must be integers for Poisson
minCNC  <- -1
maxCNC  <- 1
CNCs    <- minCNC:maxCNC # allowed outputs = gain, loss, neutral
badMapFrac <- 1/100 # allowance for some false reads mapping to a (CN+CNC)=0 bin

# initialize HMM object
message("initializing HMM object")
source(paste(createArgs$LIB_DIR, "run_HMM.R", sep="/"))
HMM <- HMM_init(
    # observation type = the properties of a genome bin, independent of an individual cell's state
    ot=list(MEDIAN=medians, # median count of bin a over all cells, scaled to the index cell's total read count
            CN=CNs),        # integer copy number of the bin for most/expected cells, from bin sizes or external training
    # observation state = the raw read count for a cell (or combination of cells) in a bin
    os=list(N=Ns), # cannot be fractional
    # hidden state = copy number change from bin baseline (NOT copy number)
    hs=list(CNC=CNCs)
)

# fill emmission probabilities
message("calculating emission probabilities")
# the output hidden states
for(CNC in CNCs){
    message(paste("   CNC: ", CNC, sep=""))
    hs_i <- HMM$hs$i[[as.character(CNC)]]
    # the observation types
    for(MEDIAN in medians){
        MEDIAN_ <- MEDIAN / 10
        for(CN in CNs){
            ot_i <- HMM$ot$i[[paste(MEDIAN, CN, sep=",")]]
            CN1N <- 1/CN * MEDIAN_ # read expected from a single allele
            exp <- min(max(MEDIAN_ + CNC * CN1N, CN1N * badMapFrac), maxN)
            # the observation states
            HMM$ep[ot_i, hs_i, Ns + 1] <- log(pmax(1e-99, dpois(Ns, exp)))
        }   
    }
}

# commit file at Rdata for future loading during segmentation
message(paste("creating model:", createArgs$HMM_FILE))
save.image(createArgs$HMM_FILE)
