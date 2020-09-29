functions{
	
	// int[] int_2D_to_1D(int[,] mat){
	// 	
	// 	int vec[num_elements(int[,1]), num_elements(int[1,])] vec;
	// 	
	// 	vec
	// 	
	// }
	matrix function_to_proportion(matrix X, matrix alpha){
		
		matrix[rows(X), cols(alpha)+1] mu;
		
		
		for(i in 1:rows(X)) {
				mu[i] = to_row_vector( softmax( to_vector( append_col([0], X[i] * alpha) ) ) );
		}
		
		return (mu);
		
	}
}
data {
	// shards
	int<lower=1> shards;
	int lv;
	
	// Reference matrix inference
	int<lower=0> GM;
	int<lower=0> S;

	// Cell types
 	int<lower=0> Q;
  int<lower=1> n_nodes;
  int<lower=1> ct_in_nodes[n_nodes];
	int<lower=1> n_levels;
  int<lower=1> ct_in_levels[n_levels];

	// reference counts
	int y[Q, GM];
	int max_y;
	matrix[ ct_in_nodes[lv], GM] ref;
	
  // Deconvolution
  int<lower=0> G_lv;
  int G_lv_linear[G_lv];

	// Lv2 tree structure parents singles
	int<lower=0> SLV2;
	int<lower=0> PLV2;
	int parents_lv2[PLV2]; // Level one parents
	int singles_lv2[SLV2]; // Level 1 leafs

	// Lv3 tree structure parents singles
	int<lower=0> SLV3;
	int<lower=0> PLV3;
	int parents_lv3[PLV3]; // Level one parents
	int singles_lv3[SLV3]; // Level 1 leafs

	// Lv4 tree structure parents singles
	int<lower=0> SLV4;
	int<lower=0> PLV4;
	int parents_lv4[PLV4]; // Level one parents
	int singles_lv4[SLV4]; // Level 1 leafs

	// Non-centered param
	real lambda_mu_prior[2];
	real lambda_sigma_prior[2];
	real lambda_skew_prior[2];
	real sigma_intercept_prior[2];

  // Proportions priors
  // lv1
  vector[ct_in_nodes[1]]  prop_1_prior[Q * (lv > 1)]; // Root

  // lv2
  vector[ct_in_nodes[2]]  prop_a_prior[Q * (lv > 2)]; // Immune cells

  // lv3
  vector[ct_in_nodes[3]]  prop_b_prior[Q * (lv > 3)]; // b cells
  vector[ct_in_nodes[4]]  prop_c_prior[Q * (lv > 3)]; // granulocyte
  vector[ct_in_nodes[5]]  prop_d_prior[Q * (lv > 3)]; // mono_derived
  vector[ct_in_nodes[6]]  prop_e_prior[Q * (lv > 3)]; // nk
  vector[ct_in_nodes[7]]  prop_f_prior[Q * (lv > 3)]; // t_cell

	// Dirichlet regression
	int A; // factors of interest
	matrix[Q,A] X;

	// Censoring
	int how_many_cens;
	int which_cens[how_many_cens];
	int which_not_cens[Q-how_many_cens];
	real<lower=0> max_unseen;
	int spt;
	real prior_survival_time[spt];

  // Local properties of the data
  vector[S] exposure_rate;
  
}
parameters {

	// Dirichlet regression
  // lv1
  matrix[A * (lv == 1) ,ct_in_nodes[1]-1]  alpha_1; // Root

	// lv2
  matrix[A * (lv == 2) ,ct_in_nodes[2]-1]  alpha_a; // Immune cells

  // lv3
  matrix[A * (lv == 3) ,ct_in_nodes[3]-1]  alpha_b; // b cells
  matrix[A * (lv == 3) ,ct_in_nodes[4]-1]  alpha_c; // granulocyte
  matrix[A * (lv == 3) ,ct_in_nodes[5]-1]  alpha_d; // mono_derived
  matrix[A * (lv == 3) ,ct_in_nodes[6]-1]  alpha_e; // natural_killer
  matrix[A * (lv == 3) ,ct_in_nodes[7]-1]  alpha_f; // t_cell

	// lv4
  matrix[A * (lv == 4) ,ct_in_nodes[8]-1]  alpha_g; // dendritic myeloid
  matrix[A * (lv == 4) ,ct_in_nodes[9]-1]  alpha_h; // macrophage
  matrix[A * (lv == 4) ,ct_in_nodes[10]-1]  alpha_i; // NK
  matrix[A * (lv == 4) ,ct_in_nodes[11]-1] alpha_l; // CD4
  matrix[A * (lv == 4) ,ct_in_nodes[12]-1] alpha_m; // CD8

	// Unknown population
	row_vector<lower=0, upper = log(max_y)>[GM] lambda_UFO;
	real<lower=0, upper=1> prop_UFO;

	// Censoring
	vector<lower=0>[how_many_cens] unseen;
	real<lower=0> prior_unseen_alpha[how_many_cens > 0];
	real prior_unseen_beta[how_many_cens > 0];

}
transformed parameters{

	matrix[Q,A] X_ = X;
	matrix[Q,A] X_scaled = X_;
	
	if(how_many_cens > 0) {
		X_[which_cens,2] = X_[which_cens,2] + unseen;
		
		// log and scale the survival days

		X_scaled[,2] = log(X_scaled[,2]);
		X_scaled[,2] = (X_scaled[,2] - mean(X_scaled[,2])) / sd(X_scaled[,2]);
	} 

}
model {
	vector[Q*GM] mu_vector;
	vector[Q*GM] sigma_vector;
	
	real sigma_intercept = 1.5;
	
	matrix[Q,  ct_in_nodes[lv]] prop = function_to_proportion(X_scaled, alpha_1);
	
	// Calculate y hat with UFO
	matrix[Q, GM] mu =
	
			// Prop matrix
			append_col(	prop * (1-prop_UFO), 	rep_vector(prop_UFO,Q) ) * 
			
			// Expression matrix
			append_row(	ref,	exp(lambda_UFO) );
	
	// Correct for exposure
	for(q in 1:Q) mu[q] = mu[q] * exposure_rate[q];
	
	// Vectorise 
	mu_vector = to_vector(mu');
	sigma_vector = 1.0 ./  exp( log(mu_vector)  * -0.4 + sigma_intercept); //((mu_vector^(-0.4)) *  exp(sigma_intercept));
	
	// Sampling
	to_array_1d(y) ~ neg_binomial_2(mu_vector, sigma_vector);

	alpha_1[1] ~ normal(0,2);
  to_vector( alpha_1[2:] ) ~ student_t(5,  0,  2.5);
  	 
	// lambda UFO
	for(i in 1:shards) lambda_UFO[i] ~ skew_normal(6.2, 3.3, -2.7);
	target += beta_lpdf(prop_UFO | 1.001, 20) * Q;

	// Censoring
	if(how_many_cens > 0){
		
		real mu_cens = prior_unseen_alpha[1] * exp(-prior_unseen_beta[1]);
		
		// unseen
		X_[which_cens,2] ~ gamma( prior_unseen_alpha[1], mu_cens);

		// Priors
		target += gamma_lpdf(X[which_not_cens,2] | prior_unseen_alpha[1], mu_cens);
	 	target += gamma_lccdf(	X[which_cens,2] | prior_unseen_alpha[1], mu_cens);
	 	
	 	// Hyperprior
	 	prior_survival_time ~ gamma( prior_unseen_alpha[1], mu_cens);
	 	
	 	prior_unseen_beta[1] ~ student_t(3, 6, 10);
  	prior_unseen_alpha[1] ~  gamma(0.01, 0.01);

	}
}
