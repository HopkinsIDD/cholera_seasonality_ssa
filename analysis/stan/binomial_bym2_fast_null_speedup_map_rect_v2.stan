// Stan model for inferring cholera occurrence seasonality 
// 
// At the lowest administrative level the model consists of a logistic regression 
// with 12 seasonal coefficients and yearly and spatial random effects. Observations
// can consists of upper level administrative units. Model extensions include
// an offset in the start of the cholera year and spatial grouping in seasonality coefficients.
// 
// This code allows for within-chain multi-threading using Stan's map_rect function as
// described here:

functions {
  
  // intFloor ----
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
  
  // interpolateLinear ----
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
  
  // count ----
  // taken from https://www.briancallander.com/posts/map_rect/speeding_up_bayesian_sampling_with_map_rect.html
  int[] count(int[] factr, int L) {
    int N = size(factr);
    int counts[L] = rep_array(0, L);
    for (i in 1:N) {
      counts[factr[i]] += 1;
    }
    return counts;
  }
  
  // prior_tau ----
  // log of the PC prior for tau (marginal spatial variance)
  real prior_tau(real tau, real alpha, real U) {
    // compute the exponential decrease paramter
    real theta = -log(alpha)/U;
    return log(theta/2) - 1.5 * log(tau) - theta * tau^(-.5);
  }
  
  // prior_rho ----
  // log of the PC prior for rho (contribution of spatial component)
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
  
  // pack_params_rect ----
  // packs the parameters into a vector to input to map_rect
  vector pack_params_rect(
    int N_lev,
    int N_units,
    int N_years,
    real[] xis,
    vector convolved_re,
    matrix etas,
    real sigma
    ) {
      int len_vec = 1 + N_lev + N_units + N_years + 1; // size of param vector
      
      vector[len_vec] param;    // param vector
      int k = 2; // counter for the parameters
      
      param[1] = len_vec;
      
      for (i in 1:N_lev) {
        param[k] = xis[i];
        k += 1;
      }
      
      for (i in 1:N_units) {
        param[k] = convolved_re[i];
        k += 1;
      }
      
      for (j in 1:N_years) {
        param[k] = etas[j, 1];
        k += 1;
      }
      
      param[k] = sigma;
      
      return(param);
    }
    
    // pack_data_rect ----
    // Packs the data into a vector for map_rect
    int[,] pack_data_rect(
      int n_offsets,
      int N_obs,
      int N_years,
      int N_years_full,
      int N_units,
      int N_lev,
      int K,
      int N_u_lower,
      int N_u_lower_yr,
      int[] starts,
      int[] ends,
      int[] map_2_u,
      int[] map_2_lev,
      int[] map_2_year,
      int[] N_obs_year, // number of observations per year
      int[] N_mapping_year, // number of mapping elements per year
      int[] uunits,
      int[] umonths,
      int[] uyears,
      int[,] offset_years,
      int[] u_years_ind,
      int[] y,
      int[] n_replicas) {
        
        int n_cnst = 8;
        int len_vec = n_cnst + N_u_lower_yr*4 + max(N_mapping_year) + max(N_obs_year)*5;    // size of data vector
        int n_shards = n_offsets * N_years;
        int year_ind [n_shards]  ;
        int offset_ind [n_shards] ;
        
        int data_array[n_shards, len_vec];
        int obs_mat[max(N_obs_year), N_years];   // extract matrix of observations per year
        int replicas_mat[max(N_obs_year), N_years]; // extract matrix of replicas per year
        int lev_mat[max(N_obs_year), N_years]; // extract matrix of replicas per year
        int starts_mat[max(N_obs_year), N_years]; // extract matrix of replicas per year
        int ends_mat[max(N_obs_year), N_years]; // extract matrix of replicas per year
        
        int mapping_mat[max(N_mapping_year), N_years]; // matrix of mapping elements
        int uunits_mat[N_u_lower_yr, N_years]; // matrix of mapping elements
        int umonths_mat[N_u_lower_yr, N_years]; // matrix of mapping elements
        int offset_year_mat[n_offsets, N_u_lower_yr, N_years];
        int k; // counter
        
        {
          int l = 1;
          for (i in 1:N_years) {
            for (j in 1:n_offsets) {
              year_ind[l] = i;
              offset_ind[l] = j;
              l += 1;
            }
          }
        }
        
        
        // extract obs and replicas matrix by year
        {
          int idx[N_years] = rep_array(1, N_years); // counters
          int map_idx[N_years] = rep_array(1, N_years); // counter of mapping indices
          for (i in 1:N_obs) {
            int yr_idx = map_2_year[i];
            obs_mat[idx[yr_idx], yr_idx] = y[i];
            replicas_mat[idx[yr_idx], yr_idx] = n_replicas[i];
            lev_mat[idx[yr_idx], yr_idx] = map_2_lev[i];
            starts_mat[idx[yr_idx], yr_idx] = map_idx[yr_idx];
            ends_mat[idx[yr_idx], yr_idx] = map_idx[yr_idx] + ends[i] - starts[i];
            
            map_idx[yr_idx] += ends[i] - starts[i] + 1;
            idx[yr_idx] += 1;
          }
        }
        
        {
          int idx[N_years] = rep_array(1, N_years); // counters
          for (i in 1:K) {
            int yr_idx = uyears[map_2_u[i]];
            mapping_mat[idx[yr_idx], yr_idx] = map_2_u[i] - (yr_idx-1) * N_units * 12; // this assumes years are consecutive in uuyears
            idx[yr_idx] += 1;
          }
        }
        
        {
          int idx[N_years] = rep_array(1, N_years); // counters
          for (i in 1:N_u_lower) {
            int yr_idx = uyears[i];
            uunits_mat[idx[yr_idx], yr_idx] = uunits[i];
            umonths_mat[idx[yr_idx], yr_idx] = umonths[i];
            for (j in 1:n_offsets) {
              offset_year_mat[j, idx[yr_idx], yr_idx] = offset_years[j, i];
            }
            idx[yr_idx] += 1;
          }
        }
        
        for (l in 1:n_shards) {
          int yr_idx = year_ind[l];
          int off_idx = offset_ind[l];
          
          k = n_cnst + 1;
          data_array[l ,1] = len_vec;
          data_array[l, 2] = N_obs_year[yr_idx];
          data_array[l, 3] = N_years;
          data_array[l, 4] = N_years_full;
          data_array[l, 5] = N_units;
          data_array[l, 6] = N_lev;
          data_array[l, 7] = N_mapping_year[yr_idx];
          data_array[l, 8] = N_u_lower_yr;
          
          for (i in 1:N_u_lower_yr) {
            data_array[l, k] = uunits_mat[i, yr_idx];
            k += 1;
          }
          
          for (i in 1:N_u_lower_yr) {
            data_array[l, k] = umonths_mat[i, yr_idx];
            k += 1;
          }
          
          for (i in 1:N_u_lower_yr) {
            data_array[l, k] = u_years_ind[yr_idx];
            k += 1;
          }
          
          for (i in 1:N_mapping_year[yr_idx]) {
            data_array[l, k] = mapping_mat[i, yr_idx];
            k += 1;
          }
          
          for (i in 1:N_obs_year[yr_idx]) {
            data_array[l, k] = lev_mat[i, yr_idx];
            k += 1;
          }
          
          for (i in 1:N_obs_year[yr_idx]) {
            data_array[l, k] = starts_mat[i, yr_idx];
            k += 1;
          }
          
          for (i in 1:N_obs_year[yr_idx]) {
            data_array[l, k] = ends_mat[i, yr_idx];
            k += 1;
          }
          
          for (i in 1:N_u_lower_yr) {
            data_array[l, k] = offset_year_mat[off_idx, i, yr_idx];
            k += 1;
          }
          
          for (i in 1:N_obs_year[yr_idx]) {
            data_array[l, k] = obs_mat[i, yr_idx];
            k += 1;
          }
          
          for (i in 1:N_obs_year[yr_idx]) {
            data_array[l, k] = replicas_mat[i, yr_idx];
            k += 1;
          }
        }
        return(data_array);
      }
      
      // compute_ll ----
      // computes the log-likelihood of a given data subset
      vector compute_ll(
        vector param,
        vector theta,
        real[] x_rs,
        int[] x_is
        ) {
          vector[1] ll; // the log likelihood to compute
          real sum_log1p;
          real sum_l;
          real x; 
          
          // data
          int len_vec = x_is[1];
          int N_obs = x_is[2];
          int N_years = x_is[3];
          int N_years_full = x_is[4];
          int N_units = x_is[5];
          int N_lev = x_is[6];
          int K = x_is[7];
          int N_u_lower = x_is[8];
          int uunits[N_u_lower];
          int umonths[N_u_lower];
          int uyears[N_u_lower];
          int map_2_u[K];
          int map_2_lev[N_obs];
          int starts[N_obs];
          int ends[N_obs];
          int off_uyears[N_u_lower];
          int y[N_obs];
          int n_replicas[N_obs];
          int k = 9;
          
          // parameters
          real xis[N_lev];
          vector[N_units] convolved_re;
          matrix[N_years_full, 1] etas;
          real sigma;
          
          // Other vars
          real ulogitp[N_u_lower];
          real prob_non_zero[N_obs];
          
          // unpack data
          for (i in 1:N_u_lower) {
            uunits[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:N_u_lower) {
            umonths[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:N_u_lower) {
            uyears[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:K) {
            map_2_u[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:N_obs) {
            map_2_lev[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:N_obs) {
            starts[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:N_obs) {
            ends[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:N_u_lower) {
            off_uyears[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:N_obs) {
            y[i] = x_is[k];
            k += 1;
          }
          
          for (i in 1:N_obs) {
            n_replicas[i] = x_is[k];
            k += 1;
          }
          
          // unpack parameters
          k = 2;
          
          for (i in 1:N_lev) {
            xis[i] = param[k];
            k += 1;
          }
          
          for (i in 1:N_units) {
            convolved_re[i] = param[k];
            k += 1;
          }
          
          for (j in 1:N_years_full) {
            etas[j, 1] = param[k];
            k += 1;
          }
          
          sigma = param[k];
          
          // compute unique lower level logitps
          for (i in 1:N_u_lower) {
            ulogitp[i] = etas[off_uyears[i], 1] + convolved_re[uunits[i]] * sigma;
          }
          
          // loop over observations and compute the probability of occurrence
          for (start_idx in 1:N_obs) {
            real xi = xis[map_2_lev[start_idx]];
            sum_log1p = 0;
            sum_l = 0;
            
            for (j in starts[start_idx]:ends[start_idx]) {
              x = xi + ulogitp[map_2_u[j]];
              sum_l += x;
              sum_log1p += log1p_exp(-1 * x);   
            }
            
            // Compute probability of occurrence
            prob_non_zero[start_idx] = log_diff_exp(sum_log1p + 1e-10, -sum_l) + sum_l;
          }
          
          // Accumulate loglikelihood
          ll[1] = binomial_logit_lpmf(y | n_replicas, prob_non_zero);
          return(ll);
        }
        
        // reduce_vec_sum ----
        // computes the sum of log-liks for each offset value across years
        vector reduce_vec_sum(
          vector ll,
          int n_offsets,
          int N_years
          ) {
            vector[n_offsets] ll_sum = rep_vector(0, n_offsets);
            int n_shards = n_offsets*N_years;
            int year_ind[n_shards] ;
            int offset_ind[n_shards] ;
            
            {
              int l = 1;
              for (i in 1:N_years) {
                for (j in 1:n_offsets) {
                  year_ind[l] = i;
                  offset_ind[l] = j;
                  l += 1;
                }
              }
            }
            
            for (i in 1:n_shards) {
              int offset_idx = offset_ind[i];
              ll_sum[offset_idx] += ll[i]; 
            }
            return(ll_sum);
          }
}

data {
  // Scalars
  int<lower=0> N_units;    // number of modeling-level admin units
  int<lower=0> N_edges;    // number of edges between admin 2 neighbours
  int<lower=1, upper=N_units> node1[N_edges];    // node1[i] adjacent to node2[i]
  int<lower=1, upper=N_units> node2[N_edges];    // and node1[i] < node2[i]
  int<lower=1> N_lev;          // number of admin levels in the data
  int<lower=0> N_obs;          // number of observations
  int<lower=0> N_years;        // number of unique years covered by the data, including offset years
  int<lower=0> N_years_obs;    // number of unique years with observations, N_years_obs <= N_years
  int<lower=0> K;              // length of observations to unique admin units x months x years mapping
  int<lower=0> N_u_lower;      // unique modelling level admin x month x year combinations (N_units x N_years x 12)
  
  // Unique combinations of admin unit x month x year
  int<lower=1, upper=12> umonths[N_u_lower];        // unique month combinations
  int<lower=1, upper=N_units> uunits[N_u_lower];    // unique admin unit combinations
  int<lower=1, upper=N_years> uyears[N_u_lower];          // unique year index combinations (1:N_years_obs)
  int<lower=1, upper=N_years> uyears_full[N_u_lower];     // unique year combinations (indices of years in u_years_ind) 
  int<lower=0, upper=N_years> u_years_ind[N_years_obs];   // unique set of indices of observed years
  
  // Mappings
  int<lower=0, upper=N_lev> map_2_lev[N_obs];       // mapping of observations to administrative level
  int<lower=0, upper=N_years> map_2_year[N_obs];    // mapping of observations to the year index¸ This setup assumes that each observation only covers a single year, so observations spanning multiple years are dropped
  int<lower=0, upper=N_u_lower> map_2_u[K];         // mapping between observation indices and admin units
  int<lower=0, upper=K> starts[N_obs];        // start index of mapping for each observation
  int<lower=0, upper=K> ends[N_obs];          // end index of mapping for each observation
  
  // Data
  int<lower=0> y[N_obs];  // count of observed cholera months at each level
  int<lower=0> n_replicas[N_obs];  // number of replications of observations for each observation
  
  // Misc for spatial re
  // The spatial component follows Riebler et al. 2016 (https://doi.org/10.1177/0962280216660421)
  // and the example of implementation in Stan in https://mc-stan.org/users/documentation/case-studies/icar_stan.html
  real<lower=0> scaling_factor;    // scales the variance of the spatial effects
  int<lower=0> N_pc_rho;           // number of elements in the PC-prior of rho
  real tabulated_rho[N_pc_rho];    // values of rho for which the prior was pre-computed 
  real tabulated_pc_prior_rho[N_pc_rho];    // pre-computed values of the prior
  real<lower=0> U;                 // PC-prior parameter U 
  real<lower=0, upper=1> alpha;    // PC-prior parameter alpha
  
  // Model specification helpers
  int<lower=0, upper =2> N_offsets;       // number of groups for which cholera year offsets are considered
}
transformed data {
  // Variables of the tabulated PC-prior for rho
  real min_rho = min(tabulated_rho);
  real max_rho = max(tabulated_rho);
  real rho_range = max_rho - min_rho;
  int n_rhos = size(tabulated_pc_prior_rho);
  real log_unif = -log(12);  // uniform prior Unif(0, 12) for the offset
  
  // Model specifications
  int<upper=N_years> offset_years[12, N_u_lower]; 
  int<lower=0, upper=1> do_offsets;   // whether to do offsets
  int<lower=0, upper=12> N_offset_elems;  // number of offset element
  
  // Additional mapping
  int<lower=0, upper=N_years> map_u_2_year[K] = uyears[map_2_u];    // mapping of the u_mapping to the year index
  
  // Additional scalars for map_rect
  int N_obs_year [N_years] = count(map_2_year, N_years);
  int N_mapping_year [N_years] = count(map_u_2_year, N_years);
  int N_u_year [N_years] = count(uyears, N_years);
  int N_u_lower_yr = N_u_year[1]; // number of unique admin x month units in each year
  
  // Variables for multi-threading with map_rect
  int n_offsets = max(N_offsets * 12, 1);    // number of offsets
  int n_shards = n_offsets * N_years_obs;    // number of data subsets for multi-threading
  real x_rs[n_shards, 0];              // empty real data vector for map_rect
  vector[0] local_params[n_shards];    // empty parameter vector for map_rect
  int n_data_rect = 8 + N_u_lower_yr*4 + max(N_mapping_year) + max(N_obs_year)*5;  // length of data vector to pass to map_rect
  int n_param_vec = 1 + N_lev + N_units + N_years + 1; // length of parameter vector to pass to map_rect
  int x_is[n_shards, n_data_rect];    // integer data vector for map_rect
  
  // Define model specifications
  if (N_offsets > 0) {
    do_offsets = 1;
    N_offset_elems = 12;
  } else {
    do_offsets = 0;
    N_offset_elems = 1;
  }
  
  // pre-define vectors of offsetted yearly random effects for each offset
  for (i in 0:11) {
    for (k in 1:N_u_lower) {
      if (do_offsets == 1) {
        offset_years[i+1, k] = (umonths[k] + i) > 12 ? uyears_full[k] + 1 : uyears_full[k];
      } else {
        offset_years[i+1, k] = uyears_full[k];
      }
    }
  }
  
  // Pack data for map_rect (function in helpers file)
  x_is = pack_data_rect(
    n_offsets,
    N_obs,
    N_years_obs,
    N_years,
    N_units,
    N_lev,
    K,
    N_u_lower,
    N_u_lower_yr,
    starts,
    ends,
    map_2_u,
    map_2_lev,
    map_2_year,
    N_obs_year,
    N_mapping_year,
    uunits,
    umonths,
    uyears,
    offset_years,
    u_years_ind,
    y,
    n_replicas);
}
parameters {
  // Yearly re
  real<lower=0> sigma_etas;  // sigma of yearly random effects
  real etas_tilde[N_years, 1];    // unscaled yearly random effects
  
  // Spatial re
  real<lower=0, upper=1> rho;    // proportion unstructured vs. spatially structured variance
  vector[N_units] theta;         // heterogeneous admin-level effects
  vector[N_units] phi;           // spatial admin-level effects
  real<lower=0> tau;         // overall marginal precision
  
  // Observation levels
  real sub_xis[N_lev-1];    // offset for the reporting of cholera at the country levels
}
transformed parameters {
  // Yearly re
  matrix[N_years, 1] etas;    // yearly random effects
  
  // Spatial re 
  vector[N_units] convolved_re = sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;   // combined scaled precision of spatial and non-spatial random effects
  real sigma = 1/sqrt(tau);        // overall marginal standard deviation
  
  // Observation levels
  real xis[N_lev] ;
  
  // For multi-threading
  vector[n_param_vec] params; // packed parameter vector for map_rect
  vector[n_shards] lls;
  vector[n_offsets] lps;
  
  xis[1] = 0;
  if (N_lev>1) {
    for (i in 2:N_lev) {
      xis[i] = sub_xis[i-1];
    }
  }
  
  // Non-centered parametrization of yearly re
  for(i in 1:N_years) {
    // scale yearly random effects
    etas[i, 1] = sigma_etas * etas_tilde[i, 1];
  }
  
  // Pack parameters for map_rect
  params = pack_params_rect(
    N_lev,
    N_units,
    N_years,
    xis,
    convolved_re,
    etas,
    sigma);
    
    // Compute log-liks with multi-threading
    lls = map_rect(compute_ll, params, local_params, x_rs, x_is);
    
    // Reduce-sum by offset
    lps = reduce_vec_sum(lls, n_offsets, N_years_obs);
}
model {
  
  // marginalize out the offset year
  target += lps;
  
  // Spatial re
  target += -0.5 * dot_self(phi[node1] - phi[node2]);
  theta ~ std_normal();
  // priors on tau and rho
  target += prior_tau(tau, alpha, U);
  target += prior_rho(logit(rho), min_rho, rho_range, tabulated_pc_prior_rho, n_rhos);
  // soft sum-to-zero constraint on phi and betas
  sum(phi) ~ normal(0, 0.001 * N_units); // equivalent to mean(phi) ~ normal(0,0.001)
  
  // Observation levels
  if (N_lev>1) {
    sub_xis ~ std_normal();
  }
  
  // Yearlry re
  sigma_etas ~ normal(0, 2.5);
  etas_tilde[, 1] ~ std_normal();
}
generated quantities {
  vector[N_obs] log_lik;    // log-likg of each observation for model comparison
  vector[N_offset_elems] ll_z1 = rep_vector(0, N_offset_elems);        // likelihood of offsets in group 1
  
  // Compute observation-level probs and log-liks for each offset
  {
    vector [N_obs] prob_non_zero[N_offset_elems, 1];  // probabilities of cholera at the upper spatial level (admin 1 or country)
    matrix[N_offset_elems, 1] lp = rep_matrix(log_unif, N_offset_elems, 1);  // likelihoods of each offset
    real sum_log1p;
    real sum_l;
    real x;
    real ulogitp[N_u_lower];
    
    for (z1 in 1:N_offset_elems) {
      int off_uyears[N_u_lower] = offset_years[z1, ]; // unique years
      
      // compute unique lower level logitps
      for (i in 1:N_u_lower) {
        ulogitp[i] =  etas[off_uyears[i], 1] + convolved_re[uunits[i]] * sigma;
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
        
        // Compute probability of occurrence
        prob_non_zero[z1, 1][start_idx] = log_diff_exp(sum_log1p + 1e-10, -sum_l) + sum_l;
      }
      
      // Accumulate loglikelihood
      lp[z1, 1] += binomial_logit_lpmf(y | n_replicas, prob_non_zero[z1, 1]);
    }
    
    // Compute offset likelihoods
    ll_z1 = softmax(to_vector(lp));
    
    // Compute observation-level likelihoods
    { 
      real log_lik_mat[N_offset_elems, 1, N_obs];
      for (i in 1:N_obs) {
        for (z1 in 1:N_offset_elems) {
          log_lik_mat[z1, 1, i] = binomial_logit_lpmf(y[i] | n_replicas[i], prob_non_zero[z1, 1][i]);
        }
        log_lik[i] = log_sum_exp(N_offsets*log_unif + to_vector(to_matrix(log_lik_mat[, , i])));
      }
    }
  }
}
