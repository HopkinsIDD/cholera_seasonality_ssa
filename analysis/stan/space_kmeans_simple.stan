data {
  int<lower=0> N;  // number of data points
  int<lower=1> D;  // number of dimensions
  int<lower=1> K;  // number of clusters
  vector[D] y[N];  // observations
  vector[D] prior_mu_location[K];
  real<lower=0> prior_mu_scale;
  int<lower=1> N_edge;
  int<lower=1,upper=N_edge> node1[N_edge];
  int<lower=1,upper=N_edge> node2[N_edge];
  real<lower=0> tau_theta; // precision of ICAR prior
  real<lower=0> sigma; // precision of ICAR prior
  
}
transformed data {
  matrix[D, D] Sigmas[K];
  
  for (k in 1:K) {
    Sigmas[k] = diag_matrix(rep_vector(sigma, D));
  }
}
parameters {
  vector[D] mus[K]; // cluster means
  simplex[K] theta[N]; // mixing parameters
}
transformed parameters {
  real<upper=0> soft_z[N, K]; // log unnormalized clusters
  
  for (n in 1:N)
  for (k in 1:K)
  soft_z[n, k] = log(theta[n, k]) + multi_normal_lpdf(y[n] | mus[k], Sigmas[k]);
}
model {
  // prior
  for (k in 1:K) {
    mus[k] ~ normal(prior_mu_location[k], prior_mu_scale);
  }
  
  for (k in 1:K) {
    // Spatial correlation
    target += -tau_theta/2 * dot_self(to_vector(logit(theta[node1, k])) - to_vector(logit(theta[node2, k])));
    // theta[n] ~ dirichlet(rep_vector(1, K));
  }
  
  // likelihood
  for (n in 1:N)
  target += log_sum_exp(soft_z[n]);
}
generated quantities {
  vector[K] probs[N];
  vector[N] log_lik;
  
  for (n in 1:N) {
    probs[n] = softmax(to_vector(soft_z[n, ]));
    log_lik[n] = log_sum_exp(soft_z[n]);
  }
}
