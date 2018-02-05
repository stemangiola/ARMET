data{
	int G;                                       // Number of marker genes
	int P;                                       // Number of cell types
	int S;                                       // Number of mix samples 
	int<lower=1> R;                              // Number of covariates (e.g., treatments)
  matrix[S,R] X;                               // Array of covariates for hierarchical regresison
  vector<lower=0, upper=1>[S] lambda[1];
	matrix<lower=0>[S,G] y;                      // Observed counts
	matrix<lower=0>[G,P] x;                      // signature counts
	vector<lower=0, upper=1>[S] p_target[1];      //This is the proportion of the whole taget -> simplex_beta * p_target
	
	real<lower=0, upper=1> theta[S];

	// Main cell types
	int<lower=1> E;
	int<lower=1> E_MU;
	int<lower=1> map_e_genes[E];
	int<lower=1> map_e_ct[E];
	int<lower=1> map_e_to_mu[E];
	vector<lower=0>[E] e_;
	int<lower=1> map_e_mu_gene[E_MU];
	int<lower=1> map_e_mu_ct[E_MU];
	int e_map_matrix_dim[P];
	real e_map_matrix_norm[P];
	int e_i_matrix[P, max(e_map_matrix_dim)];
	
	// Background
	int<lower=1> B;
	int<lower=1> B_MU;
	int<lower=1> Q ;
	int<lower=1> map_bg_genes[B] ;
	int<lower=1> map_bg_ct[B] ;
	int<lower=1> map_bg_to_mu[B];
	vector<lower=0>[B] b;
	int<lower=1> map_bg_mu_gene[B_MU];
	int<lower=1> map_bg_mu_ct[B_MU];
	matrix<lower=0, upper=1>[S,Q] beta_bg;  
	
}
transformed data{
	matrix<lower=0>[S,G] y_log;                 // log tranformation of the data
	real<lower=0> mult;
	vector<lower=1>[P] alpha_hyper_prior;

	// Estaimate gene values
	vector<lower=0>[E] e_log;                 // log tranformation of the data
	vector<lower=0>[B] b_log;                 // log tranformation of the data


	y_log = log(y+1);                           // log tranformation of the data
	mult = 2;
	if(mult<1) mult = 1;
	for(p in 1:P) alpha_hyper_prior[p] = 1;
	
	// Estaimate gene values
	e_log = log(e_+1);                           // log tranformation of the data
	b_log = log(b+1);                           // log tranformation of the data

}
parameters {
	simplex[P] beta[S];                        // Proportions,  of the signatures in the xim
	real<lower=0> sigma0[S];                       // Variance linear model 
	
	simplex[P] alpha[R];
	real<lower=0> bg_sd;
	real<lower=1> phi;
	real<lower=1> phi_phi;
	
	// Main cell types
	vector<lower=0>[E_MU] e_mu;
	vector<lower=0>[E_MU] e_sigma;
	
	// Background cell types
	vector<lower=0>[B_MU] bg_mu;
	vector<lower=0>[B_MU] bg_sigma;
	
	real<lower=0> mu_mu[2];
	real<lower=0> mu_sigma[2];
	vector<lower=0>[1] sigma_mu_e;
	real<lower=0> sigma_mu_bg;

	//vector<lower=0>[P] sigma_sigma_e;
	//real<lower=0> sigma_sigma_bg;

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

	// infer gene values
	matrix[S,G] y_hat_background;   // This is the matrix of background -> will be summed up with the matrix of the target
	matrix<lower=0>[P,G] e_mu_mat;
	matrix<lower=0>[Q,G] bg_mu_mat;

	
	for(r in 1:R) alpha_mat[r] =  logit(to_row_vector(alpha[r]));
	beta_hat =  X * alpha_mat;
	for(s in 1:S) beta_hat_hat[s] = softmax(to_vector(beta_hat[s])) * phi + 1;

	for(em in 1:E_MU) e_mu_mat[map_e_mu_ct[em], map_e_mu_gene[em]] = e_mu[em];
	for(bm in 1:B_MU) bg_mu_mat[map_bg_mu_ct[bm], map_bg_mu_gene[bm]] = bg_mu[bm];
	
	y_hat_background = beta_bg * (exp(bg_mu_mat)-1);

	for(s in 1:S) beta_target[s] = to_row_vector(beta[s]) * p_target[1,s];
	y_hat_target = beta_target * (exp(e_mu_mat)-1);
	
	// Calculate final gene value
	y_hat = y_hat_target + y_hat_background;
	
	//Mult is the multiplier inverse when the Y expression is 0
	for(s in 1:S)	sigma[s] = sigma0[s] * ((1/mult)-1) / max(log(y_hat[s]+1));   // sigma0[s] * (-mult + 1) / max(y_log[s]);
	for(s in 1:S)	sigma1[s] = sigma0[s] + sigma[s] *  max(log(y_hat[s]+1));

}
model {

	matrix[S,G] y_hat_log; 
	matrix[S,G] y_err; 
	real bern_0;
	real bern_1;


	sigma0 ~ normal(0, 0.01);
	phi_phi ~ cauchy(1,2);
	bg_sd ~ normal(0,0.01);
	phi ~ normal(0,5);
	
	// Estimation of the gene counts
	for(p in 1:P){
	int	my_map[e_map_matrix_dim[p]] =  e_i_matrix[p,1:e_map_matrix_dim[p]];
 target += normal_lpdf	(		e_log[my_map] | e_mu[map_e_to_mu[my_map]], e_sigma[map_e_to_mu[my_map]]	) * e_map_matrix_norm[p];
}
	if(min(p_target[1])<1) b_log ~ normal(bg_mu[map_bg_to_mu], bg_sigma[map_bg_to_mu]);

	e_mu ~ normal(mu_mu[1], mu_sigma[1]);
	bg_mu ~ normal(mu_mu[2], mu_sigma[2]);

	e_sigma ~ normal(sigma_mu_e[1],1);
	bg_sigma ~ normal(sigma_mu_bg,1);

	mu_mu ~ normal(5,2);
	mu_sigma ~ normal(0,2);
	sigma_mu_e ~ normal(0,1);
	sigma_mu_bg ~ normal(0,1);
	//sigma_sigma_e ~ normal(0,0.1);
	//sigma_sigma_bg ~ normal(0,0.1);

	/////////////////////////////////


	y_hat_log = log(y_hat+1);
	//for(s in 1:S) y_err[s] = sigma0[s] + y_hat_log[s] * sigma[s];
	for(s in 1:S) for(g in 1:G) y_err[s,g] = sigma0[s];

	// bern_0 = bernoulli_lpmf(0 | theta);
	// bern_1 = bernoulli_lpmf(1 | theta);
	// 
	// // Linear system
	// for(s in 1:S) for(g in 1:G){
	// 	if (y_log[s,g] == 0)
	// 		target += log_sum_exp(
	// 			bern_1, 
	// 			bern_0 + normal_lpdf(y_log[s,g] | y_hat_log[s,g],  y_err[s,g] )
	// 		);
	// 	else
	// 		target += bern_0 + normal_lpdf(y_log[s,g] | y_hat_log[s,g], y_err[s,g] );
	// }

  for(s in 1:S) y_log[s] ~ student_t( 8, y_hat_log[s],  y_err[s] ); 
  
  for(s in 1:S) beta[s] ~ dirichlet(beta_hat_hat[s]);
  for(r in 1:R) alpha[r] ~ dirichlet(alpha_hyper_prior * phi_phi);

}
generated quantities{
	vector[P] beta_gen[S];                       // Proportions,  of the signatures in the xim
	for(s in 1:S) beta_gen[s] = dirichlet_rng(beta_hat_hat[s]);
}
