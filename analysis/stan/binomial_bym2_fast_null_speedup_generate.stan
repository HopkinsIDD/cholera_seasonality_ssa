data {
  int<lower=0> N_units;    // number of modeling-level admin units
  int<lower=0> N_edges;    // number of edges between admin 2 neighbours
  int<lower=1, upper=N_units> node1[N_edges];    // node1[i] adjacent to node2[i]
  int<lower=1, upper=N_units> node2[N_edges];    // and node1[i] < node2[i]
  int<lower=1> N_lev;           // number of admin levels in the data
  // int<lower=0, upper=1> country_level; // flag whether the last level corresponds to country-level data
  int<lower=0> N_obs;    // number of observations
  int<lower=0> N_years;         // number of unique years in data
  int<lower=0> K;               // length of observations to admin nuits x months mapping
  int<lower=0> N_u_lower; // unique lower level admin unit observations (N_units x N_years x 12)
  int<lower=1, upper=N_units> uunits[N_u_lower];  // unique admin unit combinations
  int<lower=1, upper=12> umonths[N_u_lower]; // unique month combinations
  int<lower=1, upper=N_years> uyears[N_u_lower]; // unique year combinations
  int<lower=1, upper=N_years> uyears_full[N_u_lower];     // unique year combinations (indices of years in u_years_ind) 
  int<lower=0, upper=N_lev> map_2_lev[N_obs];
  int<lower=0, upper=N_u_lower> map_2_u[K];  // mapping between observation indices and admin units
  int<lower=0, upper=K> starts[N_obs];
  int<lower=0, upper=K> ends[N_obs];
  int<lower=0> y[N_obs];  // count of observed cholera months at each level
  int<lower=0> n_replicas[N_obs];  // number of replications of observations for each observation
  real<lower=0> scaling_factor; // scales the variance of the spatial effects
  int<lower=0> N_pc_rho; // number of elements in the PC-prior of rho
  real tabulated_rho[N_pc_rho];
  real tabulated_pc_prior_rho[N_pc_rho];
  real<lower=0> U;
  real<lower=0, upper=1> alpha;
  int<lower=0> N_offsets; 
}
transformed data {
  // variables of the tabulated prior for rho
  real log_unif = -log(12);  // uniform prior Unif(0, 12) for the offset
  int<upper=N_years> offset_years[12, N_u_lower]; 
  int<lower=0, upper=1> do_offsets;   // whether to do offsets
  int<lower=0, upper=12> N_offset_elems;  // number of offset element
  
  if (N_offsets > 0) {
    do_offsets = 1;
    N_offset_elems = 12;
  } else {
    do_offsets = 0;
    N_offset_elems = 1;
  }
  
  // pre-determing vectors of offsetted yearly random effects for each offset
  for (i in 0:11) {
    for (k in 1:N_u_lower) {
      if (do_offsets == 1) {
        offset_years[i+1, k] = (umonths[k] + i) > 12 ? uyears_full[k] + 1 : uyears_full[k];
      } else {
        offset_years[i+1, k] = uyears_full[k];
      }
    }
  }
}
parameters {
  real<lower=0> sigma;    // overall marginal precision
  real<lower=0> sigma_etas;  // sigma of yearly random effects
  real etas_tilde[N_years, 1] ;    // unscaled yearly random effects
  
  real sub_xis[N_lev-1];    // offset for the reporting of cholera at the country levels
  real<lower=0, upper=1> rho;    // proportion unstructured vs. spatially structured variance
  vector[N_units] theta;    // heterogeneous admin-level effects
  vector[N_units] phi;      // spatial admin-level effects
  
}
transformed parameters {
  vector[N_units] convolved_re;  // combined scaled precision of spatial and non-spatial random effects
  matrix[N_years, 1] etas;    // yearly random effects
  real ulogitp[N_u_lower];
  matrix[N_offset_elems, 1] lp = rep_matrix(log_unif, N_offset_elems, 1);  // likelihoods of each offset
  vector [N_obs] prob_non_zero[N_offset_elems, 1];  // probabilities of cholera at the upper spatial level (admin 1 or country)
  real xis[N_lev] ;
  
  xis[1] = 0;
  if (N_lev>1) {
    for (i in 2:N_lev) {
      xis[i] = sub_xis[i-1];
    }
  }
  
  for(i in 1:N_years) {
    // scale yearly random effects
    etas[i, 1] = sigma_etas * etas_tilde[i, 1];
  }
  
  // variance of each component should be approximately equal to 1
  convolved_re = sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;
  
  // Compute log-liks
  {
    real sum_log1p;
    real sum_l;
    real x; 
    
    for (z1 in 1:N_offset_elems) {
      int off_uyears[N_u_lower] = offset_years[z1, ]; // unique years
      
      // compute unique lower level logitps
      for (i in 1:N_u_lower) {
        ulogitp[i] = etas[off_uyears[i], 1] + convolved_re[uunits[i]] * sigma;
      }
      
      // loop over observations and compute the probability of occurrence
      for (start_idx in 1:N_obs) {
        real xi = xis[map_2_lev[start_idx]];
        sum_log1p = 0;
        sum_l = 0;
        
        for (k in starts[start_idx]:ends[start_idx]) {
          x = xi + ulogitp[map_2_u[k]];
          sum_l += x;
          sum_log1p += log1p_exp(-1 * x);   
        }
        prob_non_zero[z1, 1][start_idx] = log_diff_exp(sum_log1p + 1e-10, -sum_l) + sum_l;
      }
      
      // Accumulate loglikelihood
      lp[z1, 1] += binomial_logit_lpmf(y | n_replicas, prob_non_zero[z1, 1]);  
    }
  }
}
generated quantities {
  vector [N_obs] prob;    // probability of observing cholera for each observation
  vector [N_obs] sim;     // simulated observations
  vector [N_u_lower] ts_prob;    // time series of probabilities for each year and admin unit
  
  
  for (i in 1:N_u_lower) {
    ts_prob[i]  = inv_logit(ulogitp[i]);
  }
  
  for (i in 1:N_obs) {
    prob[i] = inv_logit(prob_non_zero[1,1][i]);
    
    if (prob[i] > (1 - 1E-10)) {
      prob[i] = 1 - 1E-10;
    }
    if (prob[i] < 1E-10) {
      prob[i] = 1E-10;
    }
    sim[i] = binomial_rng(n_replicas[i], prob[i]);
  }
}
