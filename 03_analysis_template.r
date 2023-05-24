# example of real data input in long format

library(rethinking)

# read in raw data

dat_raw <- read.csv("obs_work_partners_directed.csv")

#'data.frame':	5517 obs. of  8 variables:
# $ date           : chr  "2022-08-20" "2022-08-20" "2022-08-20" "2022-08-20" ...
# $ tripID         : int  12 12 12 12 12 12 12 4 4 5 ...
# $ focalID        : chr  "SB11" "SB11" "SB11" "SB11" ...
# $ senderID       : chr  "SB2" "SB11" "SB11" "SB11" ...
# $ receiverID     : chr  "SB11" "SB2" "SB5" "SB16" ...
# $ sharing        : chr  "childcare" "foraging_group_membership" "foraging_group_membership" "foraging_group_membership" ...
# $ is.baby.present: int  1 1 1 1 1 1 1 0 0 0 ...
# $ dyadID         : chr  "SB11_SB2" "SB11_SB2" "SB11_SB5" "SB11_SB16" ...

# convert ID strings to integers
IDs_send_receive <- coerce_index( dat_raw$senderID , dat_raw$receiverID )
check_index(IDs_send_receive[[1]])
check_index(IDs_send_receive[[2]])

# convert dyad ID strings to integers
dyad_ID <- coerce_index( dat_raw$dyadID )
check_index(dyad_ID)

# convert layer strings to integers
layer_ID <- coerce_index( dat_raw$sharing )
check_index(layer_ID)

# prep input data
N_obs <- nrow(dat_raw)
dat_input <- list(
    M=3, # number of layers
    N_dyads=max(dyad_ID),
    N_households=max(IDs_send_receive[[1]]),
    N_obs=N_obs,
    Y=rep(1,N_obs), # behavior is all 1 here, because zeros not recorded
    D=dyad_ID,
    HA=IDs_send_receive[[1]], # sender ID
    HB=IDs_send_receive[[2]], # receiver ID
    L=layer_ID,
    dyad_eta=2 )

m01 <- cstan( file="model_02_long.stan" , data=dat_input , chains=4 , cores=4 )

precis( m01 , depth=3 , pars=c("a","sigma_T") )
post <- extract.samples(m01)

# layer correlation matrix
round( apply(post$RHO,2:3,mean) , 2 )
#image( apply(post$RHO,2:3,mean) ,main="layer correlations")

# general give/receive matrix
round( apply(post$Rho_gr,2:3,mean) , 2 )
#image(apply(post$Rho_gr,2:3,mean),main="GR correlations")
