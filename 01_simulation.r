
# generative model of reciprocal behavior in multi-layer network

library(rethinking)

# functions to build dyadic correlation matrix for multi-layer network

# makes a nice LaTeX array
build_RHO <- function(M=2) {
    # within submatrix
    W <- diag(M)
    for ( i in 1:(M-1) ) for ( j in (i+1):M )
        W[i,j] <- concat("w_{",i,j,"}")

    # between submatrix
    R <- diag(M)
    for ( i in 1:M ) for ( j in 1:M ) {
        if ( i==j )
            R[i,j] <- concat("r_",i)
        else
            R[i,j] <- concat("r_{",concat(sort(c(i,j))),"}")
    }

    # combine
    RHO <- diag(2*M)
    RHO[1:M,1:M] <- W
    RHO[(M+1):(2*M),(M+1):(2*M)] <- W
    RHO[1:M,(M+1):(2*M)] <- R

    RHO[lower.tri(RHO)] <- ""

    return(RHO)
}

# x <- build_RHO(9)
# print(xtable(x), include.rownames = F, sanitize.text.function = identity)

# random RHO that is positive semidefinite
sim_RHO <- function(M=2,eta=5) {
    # number of parameters needed
    n <- 3*choose(M,2)
    spd <- FALSE
    while ( spd==FALSE ) {
        # sample within parameters as rlkjcorr of dim M
        w <- rlkjcorr(1,M,eta=eta)
        # sample between parameters as rlkjcorr of dim M, plus diag vector length M
        r <- rlkjcorr(1,M,eta=eta)
        diag(r) <- 2*rbeta(M,eta-1,eta-1) - 1
        # combine
        RHO <- diag(2*M)
        RHO[1:M,1:M] <- w
        RHO[(M+1):(2*M),(M+1):(2*M)] <- w
        RHO[1:M,(M+1):(2*M)] <- r
        RHO[(M+1):(2*M),1:M] <- r
        #RHO[lower.tri(RHO)] <- RHO[upper.tri(RHO)]
        require(matrixcalc)
        spd <- is.positive.semi.definite(RHO)
    }
    return(RHO)
}

# synthetic data simulation
# note for large M, eta needs to be larger or it will be hard to simulate a positive semidefinite correlation matrix for the dyad effects - the sim will hang

sim_layers <- function(
    N = 25,     # number of individuals/households
    M = 3,      # layers: number of behavior types (labor,food)
    maxt = 10,  # number of time points (days, weeks, etc)
    eta = 2,    # cholesky prior for dyad correlation sub-matrices
    Rho,        # dyad correlation matrix
    alpha,      # log base rates of each layer
    bWG,        # effect of wealth on giving in each layer
    bWR         # effect of wealth on receiving in each layer
) {

    dyads <- t(combn(N,2))
    N_dyads <- nrow(dyads)

    # now simulate directed ties for all individuals
    if ( missing(alpha) )
        alpha <- runif(M,-2,1) # base rate of ties in each domain
    y <- array(NA,c(M,N,N)) # network matrix - one network for each domain

    if ( missing(Rho) ) Rho <- sim_RHO(M,eta=eta)

    # sim tie probabilities in each dyad
    P_tie <- rmvnorm2(N_dyads, Mu=rep(0,2*M) , sigma=rep(1,2*M) , Rho=Rho )

    # sim ties
    for ( did in 1:N_dyads ) {
        A <- dyads[did,1]
        B <- dyads[did,2]
        for ( m in 1:M ) {
            # just use simulated tie strength for now, insert in adjacency matrix
            y[m,A,B] <- P_tie[ did , m ]
            y[m,B,A] <- P_tie[ did , m+M ]
        }#m
    }#i

    # simulate wealth covariate
    W <- rnorm(N) # standardized relative wealth in community
    if (missing(bWG))
        bWG <- rep(0,M) # effect of wealth on giving 
    if (missing(bWR))
        bWR <- rep(0,M) # effect of wealth on receiving 

    # now simulate behavior
    yAB <- array(NA,c(N_dyads,M,maxt))
    yBA <- array(NA,c(N_dyads,M,maxt))
    for ( i in 1:N_dyads ) {
        A <- dyads[i,1]
        B <- dyads[i,2]
        for ( j in 1:M ) {
            for ( t in 1:maxt ) {
                yAB[i,j,t] <- rpois( 1 , exp( alpha[j] + y[j,A,B] + bWG[j]*W[A] + bWR[j]*W[B] ) )
                yBA[i,j,t] <- rpois( 1 , exp( alpha[j] + y[j,B,A] + bWG[j]*W[B] + bWR[j]*W[A] ) )
            }#t time point
        }#j domain
    }#i dyad

    sim_data <- list(
        M = M,
        max_t = maxt,
        N_dyads = N_dyads,
        N_households = N,
        D = 1:N_dyads,
        HA = dyads[,1],
        HB = dyads[,2],
        YAB = yAB,
        YBA = yBA,
        dyad_eta=2,
        RHO_true=Rho,
        ties=y )

    return( sim_data )

}

# code to plot networks in each layer
if ( FALSE ) {

library(igraph)

par(mfrow=make.grid(M))

sng <- graph_from_adjacency_matrix(y[1,,])
lx <- layout_nicely(sng)
vcol <- "#DE536B"

for ( i in 1:M ) {
    sng <- graph_from_adjacency_matrix(y[i,,])
    plot( sng , layout=lx , vertex.size=8 , edge.arrow.size=0.75 , edge.width=2 , edge.curved=0.35 , vertex.color=vcol , edge.color=grau() , asp=0.9 , margin = -0.05 , vertex.label=NA  )
}

}
