
# test on synthetic data

# original adjacency matrix form

( Rho_custom <- sim_RHO( 2 , eta=2 ) )

sim_data <- sim_layers( N=25 , M=2 , eta=2 , maxt=5 , alpha=rep(0,3) , Rho=Rho_custom )

sim_dat <- sim_data
sim_dat$ties <- NULL

m01 <- cstan( file="model_02.stan" , data=sim_dat , chains=1 )

precis( m01 , depth=3 , pars=c("a","sigma_T") )
post <- extract.samples(m01)

# long form

( Rho_custom <- sim_RHO( 2 , eta=2 ) )

# sim_data <- sim_layers_long( N=25 , M=2 , eta=2 , maxt=10 , alpha=rep(0,3) , Rho=Rho_custom )

# convert sim_data to long
N_obs <- prod(dim(sim_data$YAB))*2
YY <- DD <- HA <- HB <- LL <- rep(NA,N_obs)
r <- 1
for ( i in 1:sim_data$N_dyads ) {
    for ( j in 1:sim_data$M ) {
        for ( k in 1:sim_data$max_t ) {
            YY[r] <- sim_data$YAB[i,j,k]
            DD[r] <- i
            HA[r] <- sim_data$HA[i]
            HB[r] <- sim_data$HB[i]
            LL[r] <- j
            r <- r + 1
            YY[r] <- sim_data$YBA[i,j,k]
            DD[r] <- i
            HA[r] <- sim_data$HB[i]
            HB[r] <- sim_data$HA[i]
            LL[r] <- j
            r <- r + 1
        }
    }
}

sim_dat <- list(
    M=sim_data$M,
    N_dyads=sim_data$N_dyads,
    N_households=sim_data$N_households,
    N_obs=N_obs,
    Y=YY,
    D=DD,
    HA=HA,
    HB=HB,
    L=LL,
    dyad_eta=2 )

m02 <- cstan( file="model_02_long.stan" , data=sim_dat , chains=1 )

precis( m02 , depth=3 , pars=c("a","sigma_T") )
post <- extract.samples(m02)

#
round( apply(post$RHO,2:3,mean) , 2 )
round(sim_data$RHO_true,2)
round( apply(post$Rho_gr,2:3,mean) , 2 )

par(mfrow=c(1,2))
image(sim_data$RHO_true,main="ground truth")
image(apply(post$RHO,2:3,mean),main="inferred")
