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
    int N_eta_groups,
    real[] xis,
    vector betas_g1,
    vector betas_g2,
    vector convolved_re,
    vector lambdas,
    matrix etas,
    real sigma
    ) {
      int len_vec = 1 + N_lev + N_units * 2 + 12*2 +
      N_years * N_eta_groups + 1; // size of param vector
      
      vector[len_vec] param;    // param vector
      int k = 2; // counter for the parameters
      
      param[1] = len_vec;
      
      for (i in 1:N_lev) {
        param[k] = xis[i];
        k += 1;
      }
      
      for (i in 1:12) {
        param[k] = betas_g1[i];
        k += 1;
      }
      
      for (i in 1:12) {
        param[k] = betas_g2[i];
        k += 1;
      }
      
      for (i in 1:N_units) {
        param[k] = convolved_re[i];
        k += 1;
      }
      
      for (i in 1:N_units) {
        param[k] = lambdas[i];
        k += 1;
      }
      
      for (i in 1:N_eta_groups) {
        for (j in 1:N_years) {
          param[k] = etas[j, i];
          k += 1;
        }
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
      int N_eta_groups,
      int do_grouping,
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
        
        int n_cnst = 10;
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
          data_array[l, 9] = N_eta_groups;
          data_array[l, 10] = do_grouping;
          
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
          int N_eta_groups = x_is[9];
          int do_grouping = x_is[10];
          int uunits[N_u_lower];
          int umonths[N_u_lower];
          int uyears[N_u_lower];
          int map_2_u[K];
          int map_2_lev[N_obs];
          int starts[N_obs];
          int ends[N_obs];
          int off_uyears[N_u_lower];
          int off_uyears_g2[N_u_lower];
          int y[N_obs];
          int n_replicas[N_obs];
          int k = 11;
          
          // parameters
          real xis[N_lev];
          vector[12] betas_g1;
          vector[12] betas_g2;
          vector[N_units] convolved_re;
          vector[N_units] lambdas;
          matrix[N_years_full, N_eta_groups] etas;
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
          
          if (do_grouping == 1) {
            for (i in 1:N_u_lower) { 
              off_uyears_g2[i] = off_uyears[i];
              off_uyears[i] = uyears[i];
            }
          } else {
            for (i in 1:N_u_lower) { 
              off_uyears_g2[i] = 1;
            }
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
          
          for (i in 1:12) {
            betas_g1[i] = param[k];
            k += 1;
          }
          
          for (i in 1:12) {
            betas_g2[i] = param[k];
            k += 1;
          }
          
          for (i in 1:N_units) {
            convolved_re[i] = param[k];
            k += 1;
          }
          
          for (i in 1:N_units) {
            lambdas[i] = param[k];
            k += 1;
          }
          
          for (i in 1:N_eta_groups) {
            for (j in 1:N_years_full) {
              etas[j, i] = param[k];
              k += 1;
            }
          }
          
          sigma = param[k];
          
          // compute unique lower level logitps
          for (i in 1:N_u_lower) {
            ulogitp[i] = lambdas[uunits[i]] * (betas_g1[umonths[i]] + etas[off_uyears[i], 1]) + 
            do_grouping * (1-lambdas[uunits[i]]) * (betas_g2[umonths[i]] + etas[off_uyears_g2[i], N_eta_groups]) +
            convolved_re[uunits[i]] * sigma;
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
  int<lower=1, upper=N_units> uunits[N_u_lower];    // unique admin unit combinations
  int<lower=1, upper=12> umonths[N_u_lower];        // unique month combinations
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
  int<lower=1, upper =2> N_groups;        // number of modeled seasonality groups (either 1 or 2)
  int<lower=1, upper =2> N_eta_groups;    // number of groups of yearly random effects (either 1 or 2)
  int<lower=1, upper=12> ind_betas_diff[2];    // indices where to impose ordering in seasonality coefficients
  int<lower=1, upper=12> beta_months[10];    // indices of months for which there is no ordering
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
  int<lower=0, upper=12> N_offset_elems_g2;  // number of offset elements in group 2 if grouping is done
  int<lower=0, upper=1> do_grouping;  // whether to do grouping
  int<lower=0, upper=2> N_betas_diff; // number of beta diff elements (either 0 or 2)
  int<lower=0, upper=10> N_betas_other_g2; // number of group2 beta elements other than the constrained ones (either 0 or 10)
  int<lower=0, upper=N_units> N_lambdas; 
  int<lower=1, upper=2> ind_eta_grp2;
  
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
  int n_data_rect = 10 + N_u_lower_yr*4 + max(N_mapping_year) + max(N_obs_year)*5;  // length of data vector to pass to map_rect
  int n_param_vec = 1 + N_lev + N_units * 2 + 12*2 + N_years * N_eta_groups + 1; // length of parameter vector to pass to map_rect
  int x_is[n_shards, n_data_rect];    // integer data vector for map_rect
  
  // Define model specifications
  if (N_offsets > 0) {
    do_offsets = 1;
    N_offset_elems = 12;
  } else {
    do_offsets = 0;
    N_offset_elems = 1;
  }
  
  if (N_groups > 1) {
    do_grouping = 1;
    N_betas_diff = 2;
    N_betas_other_g2 = 10;
    N_lambdas = N_units;
    
    // set indices for yearly random effects
    if (N_eta_groups > 1) {
      ind_eta_grp2 = 2;
    } else {
      ind_eta_grp2 = 1;
    }
    
  } else  {
    do_grouping = 0;
    N_betas_diff = 0;
    N_betas_other_g2 = 0;
    N_lambdas = 0;
    ind_eta_grp2 = 1;
  }
  
  // In this setup if the model has grouping then year offsets are only considered for group 2
  if(do_offsets == 1 && do_grouping == 1){
    N_offset_elems_g2 = 12;
    N_offset_elems = 1;
  } else {
    N_offset_elems_g2 = 1;
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
    N_eta_groups,
    do_grouping,
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
  real etas_tilde[N_years, N_eta_groups];    // unscaled yearly random effects
  
  // Spatial re
  real<lower=0, upper=1> rho;    // proportion unstructured vs. spatially structured variance
  vector[N_units] theta;         // heterogeneous admin-level effects
  vector[N_units] phi;           // spatial admin-level effects
  real<lower=0> tau;         // overall marginal precision
  
  // Seasonality
  vector[12] betas_g1;                // vector of monthly cholera probabilities for first group
  vector[N_betas_diff] betas_diff;    // difference between the betas between the two groups to impose ordering
  vector[N_betas_other_g2] betas_other_g2;    // other betas of group 2
  real<lower=0> tau_betas;   // precision of random walk of betas
  
  // Mixture
  vector[N_lambdas] lambdas_raw;              // spatially-correlated re for mixture paramters
  vector[N_lambdas] theta_lambdas_raw;        // spatially independent re for mixture parameters
  real<lower=0, upper=1> rho_lambdas[do_grouping];    // proportion of spatially correlated vs. uncorrelated variance in mixture paramteres
  real<lower=0> tau_lambdas[do_grouping];     // overall marginal precision of mixture parameters
  
  // Observation levels
  real sub_xis[N_lev-1];    // offset for the reporting of cholera at the country levels
}
transformed parameters {
  
  // For multi-threading
  vector[n_param_vec] params; // packed parameter vector for map_rect
  vector[n_shards] lls;
  vector[n_offsets] lps;
  
  // Yearly re
  matrix[N_years, N_eta_groups] etas;    // yearly random effects
  
  // Spatial re 
  vector[N_units] convolved_re = sqrt(1 - rho) * theta + sqrt(rho / scaling_factor) * phi;   // combined scaled precision of spatial and non-spatial random effects
  real sigma = 1/sqrt(tau);        // overall marginal standard deviation
  
  // Seasonality 
  vector[12] betas_g2 = rep_vector(1.0, 12);    // seasonality coefficients for group 2
  
  // Mixture
  vector[N_units] lambdas = rep_vector(1, N_units);    // mixture parameters
  
  // Observation levels
  real xis[N_lev] ;
  
  xis[1] = 0;
  if (N_lev>1) {
    for (i in 2:N_lev) {
      xis[i] = sub_xis[i-1];
    }
  }
  
  // Non-centered parametrization of yearly re
  for(i in 1:N_years) {
    for (l in 1:N_eta_groups) {
      // scale yearly random effects
      etas[i, l] = sigma_etas * etas_tilde[i, l];
    }
  }
  
  if (do_grouping == 1) {
    // Build spatially correlated mixture coefficients
    vector[N_units] convolved_re_lambdas;
    convolved_re_lambdas = sqrt(1 - rho_lambdas[do_grouping]) * theta_lambdas_raw + sqrt(rho_lambdas[do_grouping]/scaling_factor) * lambdas_raw;
    
    for (i in 1:N_units) {
      // mixture coefficients = invlogit of spatiall correlated and independent RE
      lambdas[i] = 1/(1+exp(-convolved_re_lambdas[i]));
    }
    
    // Build seasonality coefficients
    // constraints on ordering of betas
    betas_g2[ind_betas_diff[1]] = betas_g1[ind_betas_diff[1]] + exp(betas_diff[1]);
    betas_g2[ind_betas_diff[2]] = betas_g1[ind_betas_diff[2]] - exp(betas_diff[2]);
    
    for (i in 1:10) {
      betas_g2[beta_months[i]] = betas_other_g2[i];
    }
  }
  
  
  
  // Pack parameters for map_rect
  params = pack_params_rect(
    N_lev,
    N_units,
    N_years,
    N_eta_groups,
    xis,
    betas_g1,
    betas_g2,
    convolved_re,
    lambdas,
    etas,
    sigma);
    
    // Compute log-liks with multi-threading
    lls = map_rect(compute_ll, params, local_params, x_rs, x_is);
    
    // Reduce-sum by offset
    lps = reduce_vec_sum(lls, n_offsets, N_years_obs);
    
    
}
model {
  
  // marginalize out the offset year
  target += log_sum_exp(lps);
  
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
  
  // Seasonality
  tau_betas ~ normal(0, 2.5);
  target += 5.5 * log(tau_betas) - tau_betas/2 * (dot_self(betas_g1[2:12] - betas_g1[1:11]) + (betas_g1[1] - betas_g1[12])^2);
  sum(betas_g1) ~ normal(0, 0.001 * 12);
  
  // Yearlry re
  sigma_etas ~ normal(0, 2.5);
  for (g in 1:N_eta_groups) {
    etas_tilde[, g] ~ std_normal();
  }
  
  if (do_grouping == 1) {
    // priors on tau_lambdas
    target += prior_tau(tau_lambdas[do_grouping], alpha, U);
    // soft sum-to-zero constraint on phi and betas
    sum(lambdas_raw) ~ normal(0, 0.001 * N_units); // equivalent to mean(phi) ~ normal(0,0.001)
    theta_lambdas_raw ~ std_normal();
    // Priors for lambda_raw (up to proportionality)
    target += -0.5 * dot_self(lambdas_raw[node1] - lambdas_raw[node2]);
    // jacobian of the transformation of the ordered values of beta_g2[1] and beta_g2[7]
    target += betas_diff[1] + betas_diff[2];
    betas_diff ~ normal(log(2), .5);
    target += 5.5 * log(tau_betas) - tau_betas/2 * (dot_self(betas_g2[2:12] - betas_g2[1:11]) + (betas_g2[1] - betas_g2[12])^2);
    sum(betas_g2) ~ normal(0, 0.001 * 12);
    target += prior_rho(logit(rho_lambdas[do_grouping]), min_rho, rho_range, tabulated_pc_prior_rho, n_rhos);
  }
}
generated quantities {
  vector[N_obs] log_lik;    // log-likg of each observation for model comparison
  vector[N_offset_elems] ll_z1 = rep_vector(0, N_offset_elems);        // likelihood of offsets in group 1
  vector[N_offset_elems_g2] ll_z2 = rep_vector(0, N_offset_elems_g2);  // likelihood of offsets in group 2
  
  // Compute observation-level probs and log-liks for each offset
  {
    vector [N_obs] prob_non_zero[N_offset_elems, N_offset_elems_g2];  // probabilities of cholera at the upper spatial level (admin 1 or country)
    matrix[N_offset_elems, N_offset_elems_g2] lp = rep_matrix(log_unif, N_offset_elems, N_offset_elems_g2);  // likelihoods of each offset
    real sum_log1p;
    real sum_l;
    real x;
    real ulogitp[N_u_lower];
    
    for (z1 in 1:N_offset_elems) {
      for (z2 in 1:N_offset_elems_g2) {
        int off_uyears[N_u_lower] = offset_years[z1, ]; // unique years
        int off_uyears_g2[N_u_lower] = offset_years[z2, ]; // unique years
        
        // compute unique lower level logitps
        for (i in 1:N_u_lower) {
          ulogitp[i] = lambdas[uunits[i]] * (betas_g1[umonths[i]] + etas[off_uyears[i], 1]) +
          do_grouping * (1-lambdas[uunits[i]]) * (betas_g2[umonths[i]] + etas[off_uyears_g2[i], N_eta_groups]) +
          convolved_re[uunits[i]] * sigma;
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
          prob_non_zero[z1, z2][start_idx] = log_diff_exp(sum_log1p + 1e-10, -sum_l) + sum_l;
        }
        
        // Accumulate loglikelihood
        lp[z1, z2] += binomial_logit_lpmf(y | n_replicas, prob_non_zero[z1, z2]);
      }
    }
    
    // Compute offset likelihoods
    for (z in 1:N_offset_elems_g2) {
      ll_z1 += softmax(lp[, z]);
    }
    for (z in 1:N_offset_elems) {
      ll_z2 += softmax(to_vector(lp[z, ]));
    }
    
    // Compute observation-level likelihoods
    { 
      real log_lik_mat[N_offset_elems, N_offset_elems_g2, N_obs];
      for (i in 1:N_obs) {
        for (z1 in 1:N_offset_elems) {
          for (z2 in 1:N_offset_elems_g2) {
            log_lik_mat[z1, z2, i] = binomial_logit_lpmf(y[i] | n_replicas[i], prob_non_zero[z1, z2][i]);
          }
        }
        log_lik[i] = log_sum_exp(N_offsets*log_unif + to_vector(to_matrix(log_lik_mat[, , i])));
      }
    }
  }
}

