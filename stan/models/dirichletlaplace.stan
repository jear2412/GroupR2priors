/* /Dirichlet Laplace prior Bhattacharya 2015

Idea: replace single global scale tau by a vector of scales
(phi_1 tau, phi_2 tau, ..., phi_p tau)

Eq(9) presents a DE kernel for beta
We use the following equivalent representation since the DE might
give issues

beta_i | psi_i, phi_i, tau ~ N(0, sqrt(psi_i) phi_i*tau ) #scale parametrization
psi_i | Exp(1/2)
phi ~ Dir(alpha)
alpha= (api,...,api)
tau ~ gamma(n*api, 1/2) (shape p*api, rate 1/2)

Prop 2
api<1 gives unbounded marginals at 0


*/
data {
  int<lower=1> N; // total number of observations
  vector[N] y; // response variable
  int<lower=1> p; // number of population-level effects, includes intercept
  matrix[N, p] X; // population-level design matrix, includes a column of 1s
  // real<lower=0> sigma; // dispersion parameter
  
  //---- test data
  int<lower=1> Ntest; // total number of observations
  vector[Ntest] ytest; // test set
  matrix[Ntest, p] Xtest; // population-level design matrix including column of 1s
  
  vector<lower=0>[p - 1] alpha; // concentration vector
  real<lower=0> scale_sigma; // scale for the distribution of sigma
  
  int prior_only;  // should the likelihood be ignored?
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
  real Intercept;  // temporary intercept for centered predictors
  vector[pc] zbeta;
  //vector[pc] beta;
  simplex[pc] phi; 
  vector<lower=0>[pc] psi;
  real<lower=0> tau;
  real<lower=0> sigma;
}

transformed parameters {
  vector[pc] beta;
  vector<lower=0>[pc] lambda = phi* tau; // individual variances
  
  if(!prior_only) {
     beta= sigma*zbeta .* sqrt(lambda) .* inv_sds_X; 
  }else{
    beta= sigma*zbeta .* sqrt(lambda);
  }
  
}

model {
  
  // likelihood including constants
  if (!prior_only) {
    target += normal_id_glm_lpdf(yc | Xc, Intercept, beta, sigma);
  }
  
  // priors including constants
  target += std_normal_lpdf(zbeta); // prior for zb
 
  target += exponential_lpdf( psi | 0.5 );  //  psi ~ Exponential(1/2) 
  target += dirichlet_lpdf( phi | alpha); // phi ~ dir(alpha)
  target += gamma_lpdf(tau |  N*alpha[1], 0.5); // tau ~ gamma(api*p, 1/2) depends on N!
  target += student_t_lpdf(sigma | 3, 0, scale_sigma);
  target += normal_lpdf(Intercept | 0, 5);  // Intercept
  
}


generated quantities {
  // actual population-level intercept
  real b_Intercept = ymean+Intercept - dot_product(means_X, beta);
  
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

