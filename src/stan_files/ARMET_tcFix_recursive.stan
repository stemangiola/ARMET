data{
	int G;                                       // Number of marker genes
	int P;                                       // Number of cell types
	int S;                                       // Number of mix samples 
	int<lower=1> R;                              // Number of covariates (e.g., treatments)
  matrix[S,R] X;                               // Array of covariates for hierarchical regresison
	matrix<lower=0>[S,G] y;                      // Observed counts
	matrix<lower=0>[G,P] x;                      // signature counts
	int omit_regression;
	real<lower=0, upper=1> p_target[S];      //This is the proportion of the whole taget -> simplex_beta * p_target
	matrix[S,G] y_hat_background;   // This is the matrix of background -> will be summed up with the matrix of the target
	
	real<lower=0, upper=1> theta[S];
	int is_mix_microarray;
	real<lower=0> sigma_hyper_sd; 
	real<lower=0> phi_hyper_sd;
	real<lower=0> alpha_hyper_value;

}
transformed data{
	matrix<lower=0>[S,G] y_log;                 // log tranformation of the data
	real<lower=0> mult;
	vector<lower=1>[P] alpha_hyper;

	real bern_0;
	real bern_1;
	
	real skip_0_inflation = 0;
	real how_many_0s = 0;
	
	y_log = log(y+1);                           // log tranformation of the data
	mult = 2;
	if(mult<1) mult = 1;
	for(p in 1:P) alpha_hyper[p] = alpha_hyper_value;
	
	bern_0 = bernoulli_lpmf(0 | theta);
	bern_1 = bernoulli_lpmf(1 | theta);
	
	for(s in 1:S) for(g in 1:G) if(y[s,g]==0) how_many_0s += 1;
	if(how_many_0s/(G*S) < 0.05) skip_0_inflation = 1;
}
parameters {
	simplex[P] beta[S];                        // Proportions,  of the signatures in the xim
	real<lower=0> sigma0[S];                       // Variance linear model 
	
	simplex[P] alpha[R];
	real<lower=1> phi;
	real<lower=1> phi_phi;
	

}
transformed parameters{
	matrix[S,P] beta_target;                  // Multiply the simplex by the whole proportion of the target
	matrix[S,G] y_hat_target;                 // The predicted counts for the linear model
	matrix[S,G] y_hat; 
	real<upper=0> sigma[S];                       // Variance linear model 
	real sigma1[S];
	
	matrix[R,P] alpha_mat;                         // design matrix of the hierarchical regression
	matrix[S,P] beta_hat;                    // Predicted proportions in the hierachical linear model
	vector[P] beta_hat_hat[S];                    // Predicted proportions in the hierachical linear model

	
	for(s in 1:S) beta_target[s] = to_row_vector(beta[s]) * p_target[s];
	y_hat_target = beta_target * x';
	y_hat = y_hat_target + y_hat_background;
	
	if(is_mix_microarray==1){
		//Mult is the multiplier inverse when the Y expression is 0
		for(s in 1:S)	sigma[s] = sigma0[s] * ((1/mult)-1) / max(log(y_hat[s]+1));   // sigma0[s] * (-mult + 1) / max(y_log[s]);
		for(s in 1:S)	sigma1[s] = sigma0[s] + sigma[s] *  max(log(y_hat[s]+1));
	}
	else {
		for(s in 1:S)	sigma[s] = -0.0001;   // sigma0[s] * (-mult + 1) / max(y_log[s]);
		for(s in 1:S)	sigma1[s] = 0;
	}
	
	for(r in 1:R) alpha_mat[r] =  logit(to_row_vector(alpha[r]));
	beta_hat =  X * alpha_mat;
	for(s in 1:S) beta_hat_hat[s] = softmax(to_vector(beta_hat[s])) * phi + 1;

}
model {

	matrix[S,G] y_hat_log; 
	matrix[S,G] y_err; 

	sigma0 ~ normal(0, sigma_hyper_sd);
	phi_phi ~ cauchy(1,2);
	phi ~ normal(1,phi_hyper_sd);
	
	y_hat_log = log(y_hat+1);

	if(is_mix_microarray==1){
		for(s in 1:S) y_err[s] = sigma0[s] + y_hat_log[s] * sigma[s];
		for(s in 1:S) y_log[s] ~ normal(y_hat_log[s], y_err[s] ); 
	}
	else {
		for(s in 1:S) for(g in 1:G) y_err[s,g] = sigma0[s];

		if(skip_0_inflation == 0) 
			for(s in 1:S) for(g in 1:G){
				if (y_log[s,g] == 0)
					target += log_sum_exp(
						bern_1, 
						bern_0 + normal_lpdf(y_log[s,g] | y_hat_log[s,g],  y_err[s,g] )
					);
				else
					target += bern_0 + normal_lpdf(y_log[s,g] | y_hat_log[s,g], y_err[s,g] );
			}
		else 
			for(s in 1:S) y_log[s] ~ normal_lpdf( y_hat_log[s], y_err[s] ); 

	}
	
 
  if(omit_regression == 0) for(s in 1:S) beta[s] ~ dirichlet(beta_hat_hat[s]);
  for(r in 1:R) alpha[r] ~ dirichlet(alpha_hyper * phi_phi);

}
generated quantities{
	vector[P] beta_gen[S];                       // Proportions,  of the signatures in the xim
	for(s in 1:S) beta_gen[s] = dirichlet_rng(beta_hat_hat[s]);
}

