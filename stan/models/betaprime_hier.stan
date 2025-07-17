/* 
Hierachical Beta Prime Prior 
a ~ gamma(1, 2); 
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
  
  real<lower=0> b; // second shape of BP prior

  real<lower=0> scale_sigma; // scale for the distribution of sigma
  
  int prior_only; // should the likelihood be ignored?
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
  // local parameters 
  vector[pc] zbeta;
  real Intercept; // temporary intercept for centered predictors
  // R2D2 shrinkage parameters
  vector<lower=0, upper=1>[pc] V; // Transformation of the local variances parameter
  real<lower=0> sigma;  // dispersion parameter
  real<lower=0> a; // first shape of BP prior (recommended to be centered around 0.5)
}
transformed parameters {
  vector[pc] beta; // population-level effects
  vector<lower=0>[pc] tau2; // global scale parameter
  tau2 =  V ./ (1 - V); //scaled by sigma
  
  // compute actual regression coefficients
  beta = zbeta .* (sigma * sqrt(tau2));  
}
model {
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma); //Intercept+Xc*beta
  }
  
  // priors including constants
  //a ~ gamma(1,2);
  
  target += gamma_lpdf(a | 1, 2);
  
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
  
  // lambdas
  vector[pc] lambda = tau2; //variances of betas
  
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





