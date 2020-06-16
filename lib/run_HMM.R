
#===============================================================================
# HMM model description
#-------------------------------------------------------------------------------
# this script constructs and solves an HMM with expanded versatility

# in addition to hidden and observed _states_, it incorporates observation _types_,
# wherein emission probabilities are dynamically adjusted based on specified
# properties of observations, and a given property of an observation
# might be part of both its type and its state

# hidden and observed states and observation types can all be expanded
# to include all possible combinations of multiple discrete values

# EXAMPLE: predict the weather based on what people walking through a door are wearing
#   ot = observation types: all combinations of list(age=c('adult','child'),      observer=c('Tom','Kim'))
#   os = observed states:   all combinations of list(foowear=c('shoes','sandals'),hat=c('yes','no'))
#   hs = hidden states:     all combinations of list(temp=c('hot','cold'),        sky=c('sunny','cloudy'))
# where:
#     emission probabalities change based on who is doing the wearing AND the observing
#     age and observer have no correlation with weather, while footwear and hat do
#        e.g. children may be less likely to wear a hat in any weather (thus, age = obs type)
#        but anyone is more likely to wear a hat in cold weather       (thus, hat = obs state)
#     multiple paired measurements are made during each observation event
#     output weather states also have multiple component parts
#     in a different implementation, age might also vary with the weather
#        e.g. perhaps adults don't appear as often when it is hot (thus, age = obs type AND state)

# variable designations tend to follow:
#    http://www.cs.brown.edu/research/ai/dynamics/tutorial/Documents/HiddenMarkovModels.html
#===============================================================================


#===============================================================================
# REQUIRED: HMM model initialization
#-------------------------------------------------------------------------------
# create starting, emission and transmission probability arrays
# takes lists of values to combine for obs types and obs and hidden states
# defaults to equal starting probabilities, independent of observation type
# defaults ep and tp to zero probability, i.e. impossible, for all values
# but provides convenience shortcut to fill tp by persistence
HMM_combine_model <- function(l) {
    cmb   <- list(n_cat=length(l),
                  values=expand.grid(l, stringsAsFactors=FALSE))
    cmb$N <- nrow(cmb$values)
    cmb$i <- 1:cmb$N     
    names(cmb$i) <- apply(as.data.frame(lapply(cmb$values, as.character)), 1, paste0, collapse=",") 
    cmb
}
HMM_init <- function(ot, os, hs, p=NULL){
    HMM <- list()
    for(type in c('ot','os','hs')){
        HMM[[type]] <- HMM_combine_model(environment()[[type]])  
    }
    sp <- 1 / HMM$hs$N
    HMM$sp  <- log(array(sp, dim=c(          HMM$hs$N)))           # hidden state vector    
    HMM$ep  <- log(array(0,  dim=c(HMM$ot$N, HMM$hs$N, HMM$os$N))) # hidden vs. output state matrix, by obs type
    HMM$tp  <- log(array(0,  dim=c(          HMM$hs$N, HMM$hs$N))) # hidden vs. hidden state matrix
    if(!is.null(p)) HMM <- HMM_set_persistence(HMM, p)
    HMM # return HMM object
}
#===============================================================================


#===============================================================================
# OPTIONAL: HMM model construction, special use cases
#-------------------------------------------------------------------------------
# transition probabilities
#-------------------------------------------------------------------------------
# fill transition probabilities with a special case HMM
# where in-state transitions are preferred and all out-of-state
# transitions are equally weighted, independent of observation type,
# determined by p = persistence = probability of remaining in state
HMM_set_persistence <- function(HMM, p){
    HMM$tp[] <- log((1 - p) / (HMM$hs$N - 1)) # out-of-state transitions
    for(i in HMM$hs$i) HMM$tp[i,i] <- log(p)  # in-state transitions
    HMM
}
#-------------------------------------------------------------------------------
# NB: there is no special case for emission probabilities
# they must always be set after initialization by the caller
#-------------------------------------------------------------------------------
# set median emission probabilities over all os, stratified by ot and hs
#-------------------------------------------------------------------------------
HMM_weighted.median <- function (x, w) {
    flt <- !is.na(x) & !is.na(w)
    x   <- x[flt]
    w   <- w[flt]
    if(length(w)==0 | sum(w)==0){
        0
    } else {    
        o   <- order(x)
        x   <- x[o]
        w   <- w[o]
        p   <- cumsum(w)/sum(w)
        n   <- sum(p < 0.5)
        if (p[n + 1] > 0.5){
            x[n + 1]
        } else {
            (x[n + 1] + x[n + 2])/2
        }         
    }
}
HMM_set_median_prob <- function(HMM){ # used to estimate expected likelihoods
    HMM$wm <- log(array(0,dim=c(HMM$ot$N,HMM$hs$N,HMM$hs$N)))
    for(ot in 1:HMM$ot$N){
        for(hs_in in 1:HMM$hs$N){ # weight by the input hs
            w <- exp(HMM$ep[ot,hs_in,])
            for(hs_out in 1:HMM$hs$N){ # use ep from the output hs
                p <- exp(HMM$ep[ot,hs_out,])
                HMM$wm[ot,hs_in,hs_out] <- log(HMM_weighted.median(p, w))                
            }
        }
    }
    HMM
}
#===============================================================================


#===============================================================================
# OPTIONAL: construct observation types and states from paired observation properties
#    HMM             = model created with HMM_init
#    ot = obs_types  = list or data.frame with same value types and order as provided to HMM_init
#    os = obs_states = list or data.frame with same value types and order as provided to HMM_init
#-------------------------------------------------------------------------------
HMM_collapse_obs <- function(HMM, ot, os){
    env <- environment()
    obs <- list()    
    if(length(ot[[1]]) != length(os[[1]])){ # n obs check #1
        stop("different numbers of observations provided for ot vs. os")
    }    
    for(type in c('ot','os')){
        if(length(env[[type]]) != HMM[[type]]$n_cat){ # n obs categories check
            stop(paste("unexpected number of observation categories in", type))
        }
        if(HMM[[type]]$n_cat > 1){ # n obs check #2
            n_obs <- length(env[[type]][[1]])
            for(i in 2:HMM[[type]]$n_cat){
                if(length(env[[type]][[i]]) != n_obs){
                    stop(paste("different numbers of observations provided for", type))
                }
            }            
        }
        obs[[type]] <- HMM[[type]]$i[ if(HMM[[type]]$n_cat == 1){
            as.character(env[[type]][[1]])
        } else {
            apply(as.data.frame(lapply(env[[type]], as.character)), 1, paste0, collapse=",")
        } ]
    }   
    obs # return the observation type and state indices
}
#===============================================================================


#===============================================================================
# HMM algorithm execution
# inputs:
#    HMM  = model created with HMM_init and subsequently filled with probabilities
#    obs  = list or data.frame with ot and os columns, filled with observation indices
#           as created by HMM_convert_obs, but can also be constructed by caller
#    path = proposed set of hidden state indices for a series of observations
#-------------------------------------------------------------------------------
# common initialization of algorithms
#-------------------------------------------------------------------------------
HMM_init_algorithm <- function(HMM, env){
    env$T <- length(env$obs$ot) # length of the sequence of observations
    if(length(env$obs$os) != env$T){ # n obs check
        stop("different numbers of observations provided for ot vs. os")
    } # TODO: could put range checks on os and ot indices, but usually expect HMM_convert_obs was used
    env$is <- HMM$hs$i
    T
}
#-------------------------------------------------------------------------------
# Viterbi algorithm
#-------------------------------------------------------------------------------
# find the most likely path given a set of observations and a HMM
HMM_viterbi <- function(HMM, obs){

    # 1. initialization (observation t=1)
    HMM_init_algorithm(HMM, environment())
    delta     <- log(matrix(0,  nrow=T, ncol=HMM$hs$N))
    delta[1,] <- sapply(is, function(i) HMM$sp[i] + HMM$ep[obs$ot[[1]],i,obs$os[[1]]])    
    phi       <-     matrix(NA, nrow=T, ncol=HMM$hs$N)
    
    # 2. recursion;
    # NB: these 'for' loops are faster than apply methods with array as implemented and given recursion restrictions
    for (t in 2:T){
        pt <- t - 1
        ot_t  <- obs$ot[[t]]
        os_t  <- obs$os[[t]]
        for (j in is){     # j = this hs
            ep_j <- HMM$ep[ot_t,j,os_t]
            for (i in is){ # i = prev hs
                delta_ <- delta[pt,i] + HMM$tp[i,j] + ep_j
                if(delta[t,j] < delta_){
                    delta[t,j] <- delta_
                    phi[pt,j]  <- i
                }
            }
        }
    }
    
    # 3. termination
    prob <- -Inf
    hsi  <- rep(-1, T)
    for (j in is){
        if(prob < delta[T,j]){
            prob <- delta[T,j]
            hsi[T] <- j
        }
    }
    
    # 4. reconstruction
    for (t in (T-1):1) hsi[t] <- phi[t,hsi[t+1]]
    
    # return
    list(prob=prob, hsi=hsi) # overall likelihood and ordered set of most likely states  
}
#-------------------------------------------------------------------------------
# likelihood of specific path, using all parameters of the model
#-------------------------------------------------------------------------------
# find the likelihood of a specific path given a set of observations and a HMM
HMM_likelihood <- function(HMM, obs, hs){ # obs provided as indices, not categories
    T <- length(obs$ot) 
    if(length(obs$os) != T) stop("different numbers of observations provided for ot vs. os")
    if(length(hs)     != T) stop("different numbers of observations provided for ot vs. hs")
    if(is.character(hs)) hs <- HMM$hs$i[hs] # hs may be indices or categories
    t <- 1 # initialize first observation based on starting probabilities 
    prob <- HMM$sp[hs[t]] + HMM$ep[obs$ot[[t]],hs[t],obs$os[[t]]]
    for (t in 2:T){
        prob <- prob + HMM$tp[hs[t-1],hs[t]] + HMM$ep[obs$ot[[t]],hs[t],obs$os[[t]]]
    }
    prob
}
#-------------------------------------------------------------------------------
# observed and expected likelihoods of a specific set of hs
# using only ep, not tp
#-------------------------------------------------------------------------------
HMM_likelihood_init <- function(env){
    T <- length(env$obs$ot) 
    if(length(env$obs$os) != T) stop("different numbers of observations provided for ot vs. os")
    if(length(env$hs)     != T) stop("different numbers of observations provided for ot vs. hs")
    if(is.character(env$hs)) env$hs <- env$HMM$hs$i[env$hs] # hs may be indices or categories
    T
}
HMM_likelihood_obs <- function(HMM, obs, hs){ # likelihood IGNORING transition probabilities
    T <- HMM_likelihood_init(environment())
    sum(sapply(1:T, function(t) HMM$ep[obs$ot[[t]],hs[t],obs$os[[t]]]))
}
HMM_likelihood_exp <- function(HMM, obs, hs_in, hs){
    T <- HMM_likelihood_init(environment())
    if(is.character(hs_in)) hs_in <- HMM$hs$i[hs_in]
    lvl <- aggregate(1:T, list(ot=obs$ot, hs_in=hs_in, hs_out=hs), length)
    sum(sapply(1:nrow(lvl), function(i){
        HMM$wm[lvl[i,'ot'],lvl[i,'hs_in'],lvl[i,'hs_out']] * lvl[i,'x']        
    }))
}
#===============================================================================



##-------------------------------------------------------------------------------
## forward-backward algorithm => likelihood over model
##-------------------------------------------------------------------------------
### use the forward-backward algorithm to find the likelihood
### of a set of observations given a HMM, over all paths
##HMM_forward <- function(HMM, obs){
##
##    # 1. initialization (observation t=1)
##    T <- HMM_init_algorithm(HMM, environment())
##    
##    # 2. recursion
##    for (t in 2:T){
##        pt <- t - 1
##        ot_t  <- obs$ot[[t]]
##        os_t  <- obs$os[[t]]
##        delta[t,] <- sapply(HMM$hs$i, function(j){ # j = this hs
##            log(sum(exp(
##                sapply(HMM$hs$i, function(i) delta[pt,i] + HMM$tp[i,j]) # i = prev hs; delta is often called alpha...
##            ))) + HMM$ep[ot_t,j,os_t]
##        })
##    }
##    
##    # 3. sum to result
##    # here and above, are summing probabilities, NOT log probabilities!
##    log(sum(exp(delta[T,])))      
##}
###http://www.shokhirev.com/nikolai/abc/alg/hmm/hmm.html
###First, the probabilities for the single-symbol sequence are calculated as a
###product of initial i-th state probability and emission probability of the given
###symbol o(1) in the i-th state.
###To calculate alpha_t+1(j), we multiply
###every alpha_t(i) by the corresponding transition probability from the i-th state to
###the j-th state, sum the products over all states, and then multiply the result
###by the emission probability of the symbol o(t+1).
###Iterating the process, we can
###eventually calculate alpha_T(i), and then summing them over all states, we can obtain
###the required probability.
#
#
#HMM_likelihood_model <- function(HMM, obs){
#
#    # 1. initialization (observation t=1)
#    HMM_init_algorithm(HMM, environment())
#    alpha     <- log(matrix(0, nrow=T, ncol=HMM$hs$N))
#    alpha[1,] <- sapply(is, function(i) HMM$sp[i] + HMM$ep[obs$ot[[1]],i,obs$os[[1]]]) 
#    beta      <- log(matrix(0, nrow=T, ncol=HMM$hs$N))
#    beta[T,]  <- log(1)
#    
#    # major problem: forward-backward demands summing of probabilities, not log probabilities
#    # this will essentially always underflow on long sequences
#    # rendering this approach basically useless
#    
#    # 2. forward recursion
#    for (t in 2:T){
#        pt <- t - 1
#        ot_t  <- obs$ot[[t]]
#        os_t  <- obs$os[[t]]
#        alpha[t,] <- sapply(is, function(j){ # j = this hs
#            log(sum(exp(
#                sapply(is, function(i) alpha[pt,i] + HMM$tp[i,j]) # i = prev hs; delta is often called alpha...
#            ))) + HMM$ep[ot_t,j,os_t]
#        })
#    }
#    
#    # 3. backward recursion
#    
#    
#    # 4. 
#    
#    
#    # 3. sum to result
#    # here and above, are summing probabilities, NOT log probabilities!
#    log(sum(exp(delta[T,])))      
#}