functions {
  // linear interpolation function (intFloor and interpolateLinear) taken as is from https://github.com/stan-dev/stan/issues/1165
  int intFloor(int leftStart, int rightStart, real iReal) {
    // This is absurd. Use bisection algorithm to find int floor.
    int left;
    int right;
    left = leftStart;
    right = rightStart;
    while((left + 1) < right) {
      int mid;
      mid = left + (right - left) / 2;
      if(iReal < mid) {
        right = mid;
      } else {
        left = mid;
      }
    }
    return left;
  }
  
  // Interpolate arr using a non-integral index i
  real interpolateLinear(real[] arr, real i) {
    // Note: 1 <= i <= length(arr)
    int iLeft;
    real valLeft;
    int iRight;
    real valRight;
    // Get i, value at left. If exact time match, then return value.
    iLeft = intFloor(1, size(arr), i);
    valLeft = arr[iLeft];
    if(iLeft == i) {
      return valLeft;
    }
    // Get i, value at right.
    iRight = iLeft + 1;
    valRight = arr[iRight];
    
    // Linearly interpolate between values at left and right.
    return valLeft + (valRight - valLeft) * (i - iLeft);
  }
  
  // log of the PC prior for tau (marginal spatial variance)
  real prior_tau(real tau, real alpha, real U) {
    // compute the exponential decrease paramter
    real theta = -log(alpha)/U;
    return log(theta/2) - 1.5 * log(tau) - theta * tau^(-.5);
  }
  
  // log of the PC prior for pho (contribution of spatial component)
  real prior_rho(real rho, real min_rho, real rho_range, 
  real[] tabulated_pc_prior_rho, int nvals) {
    // prior on rho, lienar interpolation of tabulated values
    
    // fractional index based on value of rho (logit scale)
    real fraci = ((rho-min_rho)/(rho_range) * nvals) + 1.0;
    // enforce constraints on indices
    if(fraci > nvals) {fraci = nvals;}
    if(fraci < 1) {fraci = 1;}
    return interpolateLinear(tabulated_pc_prior_rho, fraci);
  }
}
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
  real min_rho = min(tabulated_rho);
  real max_rho = max(tabulated_rho);
  real rho_range = max_rho - min_rho;
  int n_rhos = size(tabulated_pc_prior_rho);
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
  real<lower=0> tau;    // overall marginal precision
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
  real<lower=0> sigma = 1/sqrt(tau);    // overall marginal precision
  
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
model {
  // marginalized likelihoods
  target += log_sum_exp(to_vector(lp));
  
  // Priors for phis and betas (up to proportionality)
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  
  if (N_lev > 1) {
    sub_xis ~ std_normal();
  }
  
  sigma_etas ~ normal(0, 2.5);
  etas_tilde[, 1] ~ std_normal();
  theta ~ std_normal();
  
  // priors on tau and rho
  target += prior_tau(tau, alpha, U);
  target += prior_rho(logit(rho), min_rho, rho_range, tabulated_pc_prior_rho, n_rhos);
  
  // soft sum-to-zero constraint on phi and betas
  sum(phi) ~ normal(0, 0.001 * N_units); // equivalent to mean(phi) ~ normal(0,0.001)
}
generated quantities {
  vector[N_obs] log_lik;
  real log_lik_mat[N_offset_elems, 1, N_obs];
  vector[N_offset_elems] ll_z1;
  
  ll_z1 = softmax(lp[, 1]);
  
  for (i in 1:N_obs) {
    for (z1 in 1:N_offset_elems) {
      log_lik_mat[z1, 1, i] = binomial_logit_lpmf(y[i] | n_replicas[i], prob_non_zero[z1, 1][i]);
    }
    log_lik[i] = log_sum_exp(N_offsets*log_unif + to_vector(to_matrix(log_lik_mat[, , i])));
  }
}

