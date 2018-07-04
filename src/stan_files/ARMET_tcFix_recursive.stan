functions{
	int do_skip_0_inflation(matrix x){
		real zeros = 0;
		int skip_0_inflation = 0;
		
		// Count the number of zeros in matrix
		for(s in 1:rows(x)) for(g in 1:cols(x)) if(y[s,g]==0) how_many_0s += 1;
		
		// If zeros are less than 5% then skip mixture model
		if(how_many_0s/(rows(x)*cols(x)) < 0.05) skip_0_inflation = 1;
		
		return skip_0_inflation;
	}
}
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
	
	real soft_prior;
	
	int P_a;
	matrix[S, P_a] p_ancestors;

}
transformed data{
	
	// Zero inflation
	real bern_0 = bernoulli_lpmf(0 | theta);
	real bern_1 = bernoulli_lpmf(1 | theta);
	int skip_0_inflation = do_skip_0_inflation(y);
	
	// Data tranformation
	matrix<lower=0>[S,G] y_log = log(y+1);                 // log tranformation of the data

}
parameters {
	simplex[P] beta[S];                        // Proportions,  of the signatures in the xim
	real<lower=0> sigma0[S];                       // Variance linear model 
	
  matrix[R,P_a - 1 + P] alpha;
	real<lower=P_a - 1 + P> phi;

}
transformed parameters{
	matrix[S,P] beta_target = (beta' * p_ancestors[P_a])';     // Multiply the simplex by the whole proportion of the target
	matrix[S,G] y_hat =  beta_target * x' + y_hat_background; 
	vector[P_a - 1 + P] beta_hat_hat[S];                        // Predicted proportions in the hierachical linear model

	for(s in 1:S) beta_hat_hat[s] = softmax( X[s] * alpha ) * phi ;
	
}
model {

	// Log transformation. Bad boy!
	matrix[S,G] y_hat_log = log(y_hat+1); 

	sigma0 ~ normal(0, sigma_hyper_sd);
	phi ~ normal(P_a - 1 + P,2);

		if(skip_0_inflation == 0) 
			for(s in 1:S) for(g in 1:G){
				if (y_log[s,g] == 0)
					target += log_sum_exp(
						bern_1, 
						bern_0 + normal_lpdf(y_log[s,g] | y_hat_log[s,g],  sigma0[s] )
					);
				else
					target += bern_0 + normal_lpdf(y_log[s,g] | y_hat_log[s,g], sigma0[s] );
			}
		else 
			for(s in 1:S) y_log[s] ~ normal_lpdf( y_hat_log[s], sigma0[s] ); 

  if(omit_regression == 0) for(s in 1:S) append_col(p_ancestors[s, 1:(P_a-1)], beta_target[s]) ~ dirichlet(beta_hat_hat[s]);
  for(r in 1:R) alpha[r] ~ normal(0,1);
  for(r in 1:R) sum( alpha[r] ) ~ normal(0,soft_prior * P);

}
generated quantities{
	vector[P_a - 1 + P] beta_gen[S];                       // Proportions,  of the signatures in the xim
	for(s in 1:S) beta_gen[s] = dirichlet_rng(beta_hat_hat[s]);
}
