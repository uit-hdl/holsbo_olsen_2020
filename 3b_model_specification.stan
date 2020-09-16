data {
  int<lower=1> N;     // number of observations
  int<lower=1> M;     // number of genes
  int<lower=1> L;     // number of detection methods
  
  int<lower=0, upper=1> metastasis[N];
  vector[N] expression[M];
  int<lower=1, upper=L> detection[N];
}

parameters {
  real mu_raw[M];
  vector[L] beta_raw[M];
  vector[L] alpha[M];
  
  real<lower=0> sd_beta_raw[M];
}

// centered parameterizations to make the sampling go better
transformed parameters {
  real<lower=0> sd_beta[M];
  real mu[M];
  vector[L] beta[M];
  
  for (m in 1:M) {
    mu[m] = .1*mu_raw[m];
    sd_beta[m] = exp(.2*sd_beta_raw[m]);
    
    for (l in 1:L) {
      beta[m][l] = mu[m] + sd_beta[m]*beta_raw[m][l];
    }
  }
}


model {
  sd_beta_raw ~ std_normal();
  mu_raw ~ std_normal();
  
  for (m in 1:M) {
    alpha[m] ~ normal(-1, 1);
    beta_raw[m] ~ std_normal();
    
    
    for (n in 1:N) {
      metastasis[n] ~ bernoulli_logit(alpha[m][detection[n]] + beta[m][detection[n]]*expression[m][n]);
    }
  }
}

generated quantities {
  // these two differ only by scaling.
  vector[L] pr_increase[M];
  vector[L] pr_increase_old[M]; // unused
  
  for (m in 1:M) {
    for (l in 1:L) {
      pr_increase_old[m][l] = inv_logit(alpha[m][l] + beta[m][l]) - inv_logit(alpha[m][l]);  // unused
      pr_increase[m][l] = inv_logit(alpha[m][l] + 0.1*beta[m][l]) - inv_logit(alpha[m][l]);  
    }
  }
}

