/* /Horseshoe prior Carvalho, Polson, Scott 2010
beta_i | lambda_i, tau \sim N(0,  sigma^2 lambda_i^2 tau^2)
lambda_i \sim HalfCauchy(0,1)
tau \sim HalfCauchy(0,sigma)

kappa_i= 1/(1+lambda_i^2)

lambda: local shrinkage
tau: global shrinkage
kappa : shrinkage factor for beta

*/
data {
  int<lower=1> N; // Number of observations
  int<lower=1> p; // Number of covariates (includes intercept)
  matrix[N, p] X; // Includes a column of 1s for intercept
  vector[N] y;
  //real<lower=0> sigma; // value for sigma
  //---- test data
  int<lower=1> Ntest; // total number of observations
  vector[Ntest] ytest; // test set
  matrix[Ntest, p] Xtest; // population-level design matrix including column of 1s
  
  real<lower=0> scale_sigma; // scale for the distribution of sigma
  
  int prior_only;  // should the likelihood be ignored?
}

transformed data {
  int pc = p - 1;
  matrix[N, pc] Xc; // centered version of X without an intercept
  vector[pc] means_X; // column means of X before centering
  vector[pc] sds_X; // sds of X before centering
  vector[pc] inv_sds_X ;
  vector[N] yc;
  real ymean;
  matrix[Ntest, pc] Xctest; // centered version of X without an intercept
  vector[pc] means_Xtest; // column means of X before centering
  
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
  real Intercept;  // temporary intercept for centered predictors
  vector[pc] zbeta;
  
  vector<lower=0>[pc] lambda;
  real<lower=0> tau;
  real<lower=0> sigma; // dispersion parameter
  
}

transformed parameters {
  vector[pc] beta;
   
  if(!prior_only) {
     beta= tau * lambda .* zbeta .* inv_sds_X; 
  } else {
     beta = tau * lambda .* zbeta;
  }
  

}

model {
  
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma); //Intercept+Xc*beta

  }
  
  // priors including constants
  target += std_normal_lpdf(zbeta);                  // zbeta ~ N(0,1)
  target += cauchy_lpdf(lambda | 0, 1);              // Half-Cauchy(0,1)
  target += cauchy_lpdf(tau | 0, sigma);             // Half-Cauchy(0, sigma)
    
  //target += std_normal_lpdf(sigma); // prior for zb
  target += student_t_lpdf(sigma | 3, 0, scale_sigma); // scale with sd(y)
            
  target += normal_lpdf(Intercept | 0, 5);  // Intercept
  
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

