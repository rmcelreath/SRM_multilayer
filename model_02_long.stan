// model for general M
// long form like the real data
data{
    int M; // number of layers
    int N_households;
    int N_dyads;
    int N_obs;
    array[N_obs] int Y; //observed behavior
    array[N_obs] int D; //dyad
    array[N_obs] int HA; //sender
    array[N_obs] int HB; //receiver
    array[N_obs] int L; //layer
    real dyad_eta; // prior concentration for dyad correlation matrix
}
parameters{
     vector[M] a;
     matrix[2*M,N_dyads] Z;
     //cholesky_factor_corr[4] L_Rho_T;
     // model the dyad effects as
     // (1) W: correlation sub-matrix for within person correlations across layers
     // (2) R: correlation sub-matrix for between person correlations across layers
     // (3) r: correlation vector of length M for between person correlations within layers
     corr_matrix[M] Rho_W;
     corr_matrix[M] Rho_R;
     vector<lower=-1,upper=1>[M] Rho_r;
     vector<lower=0>[M] sigma_T;
     // general effects
     matrix[2*M,N_households] Zgr;
     cholesky_factor_corr[2*M] L_Rho_gr;
     vector<lower=0>[2*M] sigma_gr;
}
transformed parameters{
    matrix[N_dyads,2*M] T;
    matrix[N_households,2*M] gr;
    matrix[2*M,2*M] RHO;
    matrix[M,M] W;
    matrix[M,M] R;

    // general giving/receiving effects
    gr = (diag_pre_multiply(sigma_gr, L_Rho_gr) * Zgr)';

    // dyad correlation matrix
    W = Rho_W;
    R = Rho_R;
    for ( i in 1:M ) R[i,i] = Rho_r[i];
    RHO[1:M,1:M] = W;
    RHO[(M+1):(2*M),(M+1):(2*M)] = W;
    RHO[1:M,(M+1):(2*M)] = R;
    RHO[(M+1):(2*M),1:M] = R;

    // build dyad varying effects matrix T
    // still need to mix scale parameters into this
    //L_Rho_T = cholesky_decompose(RHO);
    //T = (diag_pre_multiply( [sigma_T[1],sigma_T[2],sigma_T[1],sigma_T[2]] ,  L_Rho_T) * Z)';
    //T = ( quad_form_diag( RHO , [sigma_T[1],sigma_T[2],sigma_T[1],sigma_T[2]] ) * Z)';
    T = ( RHO * Z )';
}
model{
    vector[N_obs] lambda;
    // priors
    a ~ normal( 0 , 1 );
    // general giving/receiving
    sigma_gr ~ exponential( 1 );
    L_Rho_gr ~ lkj_corr_cholesky( 2 );
    to_vector( Zgr ) ~ normal( 0 , 1 );
    // dyads (ties)
    sigma_T ~ exponential( 1 );
    Rho_W ~ lkj_corr(dyad_eta);
    Rho_R ~ lkj_corr(dyad_eta);
    Rho_r ~ normal(0,0.5);
    to_vector( Z ) ~ normal( 0 , 1 );
    // observed variables
    for ( i in 1:N_obs ) {
        // need to fix this by knowing whether A is first or second in dyad
        // see other model for example
        // can use order of ids within dyad: smaller is always first
        if ( HA[i] < HB[i] )
            // A -> B
            lambda[i] = a[L[i]] + T[D[i], L[i]] + 
                    gr[HA[i], 1 + 2*(L[i]-1) ] + 
                    gr[HB[i], 2 + 2*(L[i]-1) ];
        else
            // B -> A (just adjust the T[] part, A/B fixed in data coding)
            lambda[i] = a[L[i]] + T[D[i], M+L[i]] + 
                    gr[HA[i], 1 + 2*(L[i]-1) ] + 
                    gr[HB[i], 2 + 2*(L[i]-1) ];
    }//i    
    Y ~ poisson_log( lambda );
}
generated quantities{
    //matrix[4,4] Rho_T;
    matrix[2*M,2*M] Rho_gr;
    Rho_gr = multiply_lower_tri_self_transpose(L_Rho_gr);
    //Rho_T = multiply_lower_tri_self_transpose(L_Rho_T);
}
