// Bayesian Time-Varying Site-Occupancy Model with Probit Link
// Seth Temple, sdtemple@lanl.gov
// August 31, 2020

data {
    int<lower=0> nsites;
    int<lower=0> nseasons;
    int<lower=0> nperiods;
    int<lower=0> yna[nsites, nseasons, nperiods,1];
    int<lower=0, upper=1> na[nsites, nseasons, nperiods,1]; // indicator array for missing values
    int<lower=1> nbeta; // # of site-occupancy effects
    row_vector[nbeta] X[nsites, nseasons, nperiods];
    int<lower=1> nalpha; // # of detection effects
    row_vector[nalpha] W[nsites, nseasons, nperiods];
}

parameters {
    vector[nbeta] beta;
    vector[nalpha] alpha;
}

transformed parameters {
    // declarations
    real<lower=0, upper=1> psi[nsites, nseasons, nperiods];
    real<lower=0, upper=1> p[nsites, nseasons, nperiods];

    // define probit relationship
    for(i in 1:nsites){
        for(j in 1:nseasons){
            for(k in 1:nperiods){
                psi[i,j,k] = Phi(X[i,j,k] * beta);
            }
        }
    }
    
    for(i in 1:nsites){
        for(j in 1:nseasons){
            for(k in 1:nperiods){
                p[i,j,k] = Phi(W[i,j,k] * alpha);
            }
        }
    }
}

model {
    // priors
    beta ~ normal(0,1);
    alpha ~ normal(0,1);

    // likelihood
    for(i in 1:nsites){
        for(j in 1:nseasons){
            for(k in 1:nperiods){
                if(na[i,j,k,1] == 1){
                    if(yna[i,j,k,1] == 1){
                        target += log(psi[i,j,k]);
                        target += bernoulli_lpmf(yna[i,j,k,1] | p[i,j,k]);
                    } else{
                        target += log_sum_exp(log(psi[i,j,k]) + bernoulli_lpmf(yna[i,j,k,1] | p[i,j,k]), log(1 - psi[i,j,k]));
                    }
                }
            }  
        }
    }
}