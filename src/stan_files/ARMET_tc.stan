data {

	// Reference matrix inference

	int<lower=0> G;
	int<lower=0> S;
	int CL;

 	int counts_linear[CL] ;
	int G_linear[CL] ;
	int S_linear[CL] ;

	real<upper=0> sigma_slope;
  real<lower=0>sigma_sigma;


}
parameters {

  // Gene-wise properties of the data
  vector[G] lambda_log;

  vector[S] exposure_rate;

	real<lower=0> lambda_mu;
  real<lower=0> lambda_sigma;
  real<upper=0> lambda_skew;

  real sigma_intercept;


}

model {

  // Overall properties of the data
  lambda_mu ~ normal(0,2);
	lambda_sigma ~ normal(0,2);
	lambda_skew ~ normal(0,1);
	sigma_intercept ~ student_t(8, 0, 1);

	exposure_rate ~ normal(0,1);
	sum(exposure_rate) ~ normal(0, 0.001 * S);
	lambda_log ~ skew_normal(lambda_mu, lambda_sigma, lambda_skew);

	counts_linear ~ neg_binomial_2_log(lambda_log[G_linear] + exposure_rate[S_linear], 1.0 ./ exp( sigma_slope * lambda_log[G_linear] + sigma_intercept));
	//counts_linear ~ poisson_log(lambda_log[G_linear] );

}
