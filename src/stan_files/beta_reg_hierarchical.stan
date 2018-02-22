data{
	int I;																			// Number of iterations in previous model
	int K;                                      // Number of groups
	int N;                                      // Number of observations
	int<lower=1> R;                             // Number of covariates 
	matrix[N,R] X;                              // Design matrix
	vector[I] beta_posterior[N, K];                       // Proportions posterior

}
parameters {

	vector<lower=0, upper=1>[K] beta_mu[N];                       	// Proportions
	vector<lower=0>[K] beta_sd[N];                       	// Proportions
	real<lower=0.0001> phi2; 
  matrix[R,K] alpha2; 
	real<lower=0> bg_sd2[R];                       	// Variance of the prior distribution to alpha

}
transformed parameters{
	matrix[N, K] beta_hat2; 
	beta_hat2 = inv_logit(X * alpha2); 
}
model {
	
	// Priors
	phi2 ~ cauchy(0,2.5);
	bg_sd2 ~ normal(0,1);

	// beta 
	for(n in 1:N) for(k in 1:K) beta_posterior[n,k] ~ normal(beta_mu[n,k], beta_sd[n,k]);
	
	// Linear model for beta regression
	for(n in 1:N) beta_mu[n] ~ beta(beta_hat2[n] * phi2, (1 - beta_hat2[n]) * phi2); 
	
	// Prior distribution on the background cluster
	for(r in 1:R) alpha2[r] ~ normal(0, bg_sd2[r]);
	
}


