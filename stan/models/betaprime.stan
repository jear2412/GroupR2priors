/* 
Beta Prime Prior
Bai 2021 
tau^2 ~ BetaPrime(a, b)
b_i ~ N(0, sigma^2 tau_i^2)

*/
data {
  int<lower=1> N; // total number of observations
  vector[N] y; // response variable
  int<lower=1> p; // number of population-level effects, includes intercept
  matrix[N, p] X; // population-level design matrix, includes a column of 1s
  //real<lower=0> sigma; // dispersion parameter
  
  //---- test data
  int<lower=1> Ntest; // total number of observations
  vector[Ntest] ytest; // test set
  matrix[Ntest, p] Xtest; // population-level design matrix including column of 1s
  
  // data for the BP prior
  real<lower=0> a; // first shape of BP prior (recommended 0.5 for shrinkage)
  real<lower=0> b; // second shape of BP prior
  
  real<lower=0> scale_sigma; // scale for the distribution of sigma
  
  int prior_only; // should the likelihood be ignored?
}
transformed data {
  int pc = p - 1;
  matrix[N, pc] Xc; // centered version of X without an intercept
  vector[pc] means_X; // column means of X before centering
  vector[pc] sds_X; // sds of X before centering
  vector[pc] inv_sds_X; 
  
  vector[N] yc;
  real ymean;
  
  for (i in 2 : p) {
    means_X[i - 1] = mean(X[ : , i]);
    sds_X[i - 1] = sd(X[ : , i]);
    inv_sds_X[i-1] = 1.0 / sds_X[i-1];
    Xc[ : , i - 1] = X[ : , i] - means_X[i - 1];
  }
  
  
  ymean = mean(y);
  for (i in 1 : N) {
    yc[i] = y[i] - ymean;
  }
  
}
parameters {
  // local parameters 
  vector[pc] zbeta;
  real Intercept; // temporary intercept for centered predictors
  // BP parameters
  vector<lower=0, upper=1>[pc] V; // Transformation of the local variance
  real<lower=0> sigma;  // dispersion parameter
}
transformed parameters {
  vector[pc] beta; // population-level effects
  vector<lower=0>[pc] lambda2; // local variances parameter
  lambda2 = V ./ (1 - V); //scaled by sigma
  
  // compute actual regression coefficients
  if(!prior_only){
     beta = zbeta .* (sigma * sqrt(lambda2)) .* inv_sds_X  ;  
  }else{
     beta = zbeta .* (sigma * sqrt(lambda2)) ;  
  }
  
}

model {
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma); //Intercept+Xc*beta
  }
  
  // priors including constants
  
  //V ~ Beta(a,b)
  target += beta_lpdf(V | a, b);
  
  target += std_normal_lpdf(zbeta); //normal distribution zbeta
  
  target += student_t_lpdf(sigma | 3, 0, scale_sigma ); // notice that the t distribution is scaled
  
  target += normal_lpdf(Intercept | 0, 5); // Intercept  
}
generated quantities {
  // actual population-level intercept
  real b_Intercept = ymean + Intercept - dot_product(means_X, beta);
  
  vector[N] log_lik;
  array[N] real y_tilde;
  vector[N] mu_tilde = rep_vector(0.0, N) +  b_Intercept + X[, 2:p] * beta;
  
  vector[Ntest] log_lik_test;
  array[Ntest] real y_tilde_test;
  vector[Ntest] mu_tilde_test = rep_vector(0.0, Ntest) + b_Intercept + Xtest[,2:p] * beta;

  //Predictive R2
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





