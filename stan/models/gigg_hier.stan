// GIGG prior

// ag ~ Ga(1,2)
// bg ~ Ga(1,2)

// gamma_g ~ Ga(ag, 1)
// lambda_gj ~ Inverse Gamma(gg, 1)
// gamma_g*lambda_gj ~ BetaPrime(ag, bg)
// tau2 ~ student_t(0, sigma)
// sigma ~ student_t(3,0, scale_sigma)

// beta ~ N(0, gamma_g*lambda_gj)

data {
  int<lower=1> N; // number of observations
  vector[N] y; // observations
  int<lower=1> p; // total number of covariates (includes intercept)
  int<lower=0> G; // number of groups
  array[p - 1] int pg; // vector of group sizes (of length p-1) // alternative array[p] int pg;
  matrix[N, p] X; //total covariate matrix (includes column of 1s)
  
  //real<lower=0> sigma; // dispersion parameter
  
  //---- test data
  int<lower=1> Ntest; // total number of observations
  vector[Ntest] ytest; // test set
  matrix[Ntest, p] Xtest; // population-level design matrix including column of 1s

  real<lower=0> scale_sigma; // scale for the distribution of sigma
  vector<lower=0>[G] bg;
  int prior_only;  // should the likelihood be ignored?
}
transformed data {
  int pc = p - 1;
  matrix[N, pc] Xc; // centered version of X without an intercept
  vector[pc] means_X; // column means of X before centering
  vector[pc] sds_X; // sds of X before centering
  vector[N] yc;
  real ymean;
  matrix[Ntest, pc] Xctest; // centered version of X without an intercept
  vector[pc] means_Xtest; // column means of X before centering
  
  for (i in 2 : p) {
    means_X[i - 1] = mean(X[ : , i]);
    sds_X[i - 1] = sd(X[ : , i]);
    //Xc[, i - 1] = (X[, i] - means_X[i - 1])/sds_X[i - 1];
    Xc[ : , i - 1] = X[ : , i] - means_X[i - 1];
  }
  
  ymean = mean(y);
  for (i in 1 : N) {
    yc[i] = y[i] - mean(y);
  }
  
  for (i in 2 : p) {
    means_Xtest[i - 1] = mean(Xtest[ : , i]);
    Xctest[ : , i - 1] = Xtest[ : , i] - means_Xtest[i - 1];
  }

}
parameters {
  //real alpha; // intercept
  real Intercept; // temporary intercept for centered predictors
  real<lower=0> sigma; // noise
  
  // Parameters 
  real<lower=0> tau2; // global shrinkage parameter
  vector<lower=0>[G] gamma2; // group shrinkage factor
  vector<lower=0>[pc] lambda2; // local shrinkage parameter
  vector<lower=0>[G] ag;
  //vector<lower=0>[G] bg;
  
  // Group level terms
  vector[pc] z_beta;
}
transformed parameters {
  vector<lower=0>[pc] Sigma; // vector of variances
  vector[pc] beta;
  
  for (j in 1 : pc) {
    Sigma[j] = tau2 * gamma2[pg[j]] * lambda2[j];
  }

  
  beta= z_beta .* sqrt(Sigma); 
  
  
}
model {
  // priors
  Intercept ~ normal(0, 5);
  z_beta ~ std_normal();
  tau2 ~ student_t(1, 0, sigma); // TODO: Should we change this? Discuss
  sigma ~ student_t(3, 0, scale_sigma);

  gamma2 ~ gamma(ag, 1); 
  
  for (j in 1:pc){
  lambda2[j] ~ inv_gamma(bg[pg[j]], 1);
  }
  
  ag ~ gamma(1,2);
  //bg ~ gamma(1,2);
  
  // likelihood
  if (!prior_only) {
    yc ~ normal_id_glm(Xc, Intercept, beta, sigma);
  }
  
  
}
generated quantities {
  real b_Intercept = ymean+Intercept - dot_product(means_X, beta);

  
  vector[N] log_lik;
  array[N] real y_tilde;
  vector[N] mu_tilde = rep_vector(0.0, N) +  b_Intercept + X[, 2:p] * beta;
  
  vector[Ntest] log_lik_test;
  array[Ntest] real y_tilde_test;
  vector[Ntest] mu_tilde_test = rep_vector(0.0, Ntest) + b_Intercept + Xtest[,2:p] * beta;
  
  //--- R2
  real<lower=0, upper=1> pred_R2 = variance(mu_tilde) / (variance(mu_tilde) + sigma^2 );
  real<lower=0, upper=1> pred_R2_test = variance(mu_tilde_test) / (variance(mu_tilde_test) + sigma^2 );
  
  //---y_tilde calc
  for (n in 1 : N) {
    log_lik[n] = normal_lpdf(y[n] | mu_tilde[n], sigma);
    y_tilde[n] = normal_rng(mu_tilde[n], sigma); //copy and paste model (executed once per sample) 
  }
  
  //---y_tilde test calc
  for (n in 1 : Ntest) {
    log_lik_test[n] = normal_lpdf(ytest[n] | mu_tilde_test[n], sigma);
    y_tilde_test[n] = normal_rng(mu_tilde_test[n], sigma); //copy and paste model (executed once per sample) 
  }
}


