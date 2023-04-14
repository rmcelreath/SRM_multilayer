
# test on synthetic data

( Rho_custom <- sim_RHO( 3 , eta=2 )
sim_data <- sim_layers( N=25 , M=3 , eta=2 , alpha=rep(0,3) , Rho=Rho_custom )

sim_dat <- sim_data
sim_dat$ties <- NULL

m01 <- cstan( file="model_02.stan" , data=sim_dat , chains=1 )

precis( m01 , depth=3 , pars=c("a","sigma_T") )

post <- extract.samples(m01)
round( apply(post$RHO,2:3,mean) , 2 )
round(sim_data$RHO_true,2)
round( apply(post$Rho_gr,2:3,mean) , 2 )

par(mfrow=c(1,2))
image(sim_data$RHO_true,main="ground truth")
image(apply(post$RHO,2:3,mean),main="inferred")
