functions{
	int do_skip_0_inflation(matrix x){
		real zeros = 0;
		int skip_0_inflation = 0;
		int how_many_0s = 0;

		// Count the number of zeros in matrix
		for(s in 1:rows(x)) for(g in 1:cols(x)) if(x[s,g]==0) how_many_0s += 1;

		// If zeros are less than 5% then skip mixture model
		if(how_many_0s * 1.0 / (rows(x)*cols(x)) < 0.05) skip_0_inflation = 1;

		return skip_0_inflation;
	}

	matrix vector_array_to_matrix(vector[] x) {
		matrix[size(x), rows(x[1])] y;
		for (m in 1:size(x))
		  y[m] = x[m]';
		return y;
	}

	 vector reg_horseshoe(
					vector zb,
					real aux1_global ,
					real aux2_global,
					vector  aux1_local ,
					vector aux2_local ,
					real  caux,
					real scale_global ,
					real slab_scale
					) {
    int K = rows(zb);

    // Horseshoe variables
		real tau ; // global shrinkage parameter
		vector [ K] lambda ; // local shrinkage parameter
		vector [ K] lambda_tilde ; // ’ truncated ’ local shrinkage parameter
		real c; // slab scale

		// Horseshoe calculation
		lambda = aux1_local .* sqrt ( aux2_local );
		tau = aux1_global * sqrt ( aux2_global ) * scale_global * 1 ;
		c = slab_scale * sqrt ( caux );
		lambda_tilde = sqrt ( c ^2 * square ( lambda ) ./ (c ^2 + tau ^2* square ( lambda )) );
		return  zb .* lambda_tilde * tau ;
  }

}
data{

	// Deconvolution
	int G;                                       // Number of marker genes
	int P;                                       // Number of cell types
	int S;                                       // Number of mix samples
	int is_mix_microarray;
	matrix<lower=0>[S,G] y;                      // Observed counts
	matrix<lower=0>[G,P] x;                      // signature counts
	matrix[S,G] y_hat_background;                // This is the matrix of background
	real<lower=0, upper=1> theta[S];

	// Regression
	int<lower=1> R;                              // Number of covariates (e.g., treatments)
	int P_a;
  matrix[S,R] X;                               // Array of covariates for hierarchical regresison
	int omit_regression;
	matrix[S, P_a] p_ancestors;

	// Horseshoe
	real < lower =0 > par_ratio ; // proportion of 0s
	real < lower =1 > nu_global ; // degrees of freedom for the half -t prior
	real < lower =1 > nu_local ; // degrees of freedom for the half - t priors
	real < lower =0 > slab_scale ; // slab scale for the regularized horseshoe
	real < lower =0 > slab_df; // slab degrees of freedom for the regularized

}
transformed data{

	// Transformed dimensions
	int P_tot = P_a - 1 + P;

	// Zero inflation
	real bern_0 = bernoulli_lpmf(0 | theta);
	real bern_1 = bernoulli_lpmf(1 | theta);
	int skip_0_inflation = do_skip_0_inflation(y);

	// Data log tranformation
	matrix<lower=0>[S,G] y_log = log(y+1);

	// Horseshoe
	real < lower =0 > scale_global = 1; //par_ratio / sqrt(1.0 * S); // scale for the half -t prior for tau

	vector<lower=0>[S] sigma = rep_vector(0.5, S);

print(skip_0_inflation);

}
parameters {

	// Deconvolution
	simplex[P] beta[S];                // Proportions,  of the signatures in the xim
	vector<lower=0, upper=1>[S] beta_dump;
	//vector<lower=0>[S] sigma;           // Variance linear model
	matrix<lower=0>[G,1] x_dump;
	// Regression
	matrix[R,P_tot] extrinsic_raw;
	real<lower=0> phi_raw;

	// Horseshoe
	real < lower =0 > aux1_global ;
	real < lower =0 > aux2_global ;
	vector < lower =0 >[P_tot] aux1_local ;
	vector < lower =0 >[P_tot] aux2_local ;
	real < lower =0 > caux ;
}
transformed parameters{

	real<lower=0> phi = inv(sqrt(phi_raw));  	// Variance Dirichlet
	matrix[R,P_tot] extrinsic;             	// Covariates
	matrix[S,P_tot] beta_global;
	vector[P_tot] beta_hat_hat[S];    // Predicted proportions in the hierachical linear model


	matrix[S,G] y_hat =
		// Beta target
		(
			// Insert dump cell type
			append_row(
				vector_array_to_matrix(beta)' *
				diag_matrix(1-beta_dump),
				to_row_vector(beta_dump)
			) *
			// Multiplication by row for target
			diag_matrix(p_ancestors[,P_a])
		)' *
		// Append the dump ct gene values
		append_col(x, x_dump)' +
		// Sum the background gene values
		y_hat_background;

	// Building matrix factors of interest
	extrinsic = extrinsic_raw;
	if(R > 1) extrinsic[2] = to_row_vector(
		reg_horseshoe(
			to_vector(extrinsic_raw[2]),
			aux1_global ,
			aux2_global,
			aux1_local ,
			aux2_local ,
			caux,
			scale_global,
			slab_scale
		)
	);

	// Regression
	for(s in 1:S) beta_hat_hat[s] = softmax( to_vector( X[s] * extrinsic ) ) * phi + 1;
	beta_global =
		append_col(
			// Matrix of ancestor proportions
			p_ancestors[, 1:(P_a-1)],
			// Beta target
			(
				vector_array_to_matrix(beta)' *
				diag_matrix(p_ancestors[,P_a])
			)'
		);

}
model {

	// Log transformation. Bad boy!
	matrix[S,G] y_hat_log = log(y_hat+1);

// Horseshoe
	aux1_local ~ normal (0 , 1);
	aux2_local ~ inv_gamma (0.5* nu_local , 0.5* nu_local );
	aux1_global ~ normal (0 , 1);
	aux2_global ~ inv_gamma (0.5* nu_global , 0.5* nu_global );
	caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df );

	// Deconvollution
	//sigma ~ normal(0.5, 0.1);

	// If 0 inflation
	if(skip_0_inflation == 0)
		for(s in 1:S) for(g in 1:G){
			if (y_log[s,g] == 0)
				target += log_sum_exp(
					bern_1,
					bern_0 + student_t_lpdf(y_log[s,g] |10, y_hat_log[s,g],  sigma[s] )
				);
			else
				target += bern_0 + student_t_lpdf( y_log[s,g] |10, y_hat_log[s,g], sigma[s] );
		}

	// If not 0 inflation
	else for(s in 1:S) y_log[s] ~ student_t(10, y_hat_log[s], sigma[s] );

	// Regression
	x_dump[1:G,1] ~ lognormal(0,1);
	beta_dump ~ beta(2, 10);
	phi_raw ~ normal(0,1);
	extrinsic_raw[1] ~ normal(0,5);
	sum(extrinsic_raw[1]) ~ normal(0,0.01 * P_tot) ;
	if(R > 1) extrinsic_raw[2] ~ normal (0 , 1);
	if(R > 2) for(r in 3:R) extrinsic_raw[r]  ~ normal(0,1);
  if(omit_regression == 0)
  	for(s in 1:S) to_vector( beta_global[s] ) ~ dirichlet(beta_hat_hat[s]);

}
generated quantities{
	vector[P_tot] beta_gen[S];                       // Proportions,  of the signatures in the xim
	for(s in 1:S) beta_gen[s] = dirichlet_rng(beta_hat_hat[s]);
}
