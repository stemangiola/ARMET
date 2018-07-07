data{
	int S;
	int P;                                      // Number of groups
	int N;                                      // Number of observations
	int<lower=1> R;                             // Number of covariates
	matrix[S,R] X;                              // Design matrix
	vector<lower=0, upper=1>[N] beta[S,P];
}
parameters {
	matrix<lower=0, upper=1>[S,P] beta_a;
	matrix<lower=0>[S,P] beta_b;
	matrix<lower=0, upper=1>[S,P] beta_hat;
	//
	real<lower=0> phi;
	matrix[R,P] intrinsic;
	real<lower=0> sigma;                       	// Variance of the prior distribution to alpha

}
transformed parameters{
	matrix[S, P] beta_mu = inv_logit(X * intrinsic);
	matrix<lower=0>[S,P] beta_b_inv = inv(beta_b);
}
model {

	// Process posterior of previous model
	for(s in 1:S) for(p in 1:P) beta[s,p] ~ beta(beta_b_inv[s,p] * beta_a[s,p], beta_b_inv[s,p] * (1 - beta_a[s,p]));
	for(s in 1:S) beta_a[s] ~ beta(5,5);
	for(s in 1:S) beta_b[s] ~ normal(0,1);
	for(s in 1:S) for(p in 1:P) beta_hat[s,p] ~ beta(beta_b_inv[s,p] * beta_a[s,p], beta_b_inv[s,p] * (1 - beta_a[s,p]));

	// Likelihood
	for(s in 1:S) beta_hat[s] ~ beta(beta_mu[s] * phi, (1 - beta_mu[s]) * phi);

	// Priors
	phi ~ cauchy(0,2.5);
	sigma ~ cauchy(0,2.5);

	// Prior distribution on the background cluster
	intrinsic[1] ~ normal(1.0/P, sigma);
	intrinsic[2] ~ normal(0, 2);
	if(R > 2) for(r in 3:R) intrinsic[r]  ~ normal(0,1);

}
generated quantities{
	matrix[S,P] beta_gen;
	for(s in 1:S) for(p in 1:P) beta_gen[s,p] = beta_rng(beta_b_inv[s,p] * beta_a[s,p], beta_b_inv[s,p] * (1 - beta_a[s,p]));
}



