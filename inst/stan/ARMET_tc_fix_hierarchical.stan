functions{

matrix vector_array_to_matrix(vector[] x) {
		matrix[size(x), rows(x[1])] y;
		for (m in 1:size(x))
		  y[m] = x[m]';
		return y;
}

vector[] which(int x, vector[] a, vector[] b, vector[] c, vector[] d){
		if(x == 1) return(a);
		if(x == 2) return(b);
		if(x == 3) return(c);
		else return(d);
	}
	
vector[] append_vector_array(vector[] v1, vector[] v2){
	vector[num_elements(v1[1]) + num_elements(v2[1])] v3[num_elements(v1[,1])];

	v3[,1:num_elements(v1[1])] = v1;
	v3[,num_elements(v1[1])+1:num_elements(v3[1])] = v2;

	return v3;
}

vector[] multiply_by_column(vector[] v, real[] r){
	int n_rows = num_elements(v[,1]);
	int n_cols = num_elements(v[1]);

	vector[n_cols] v_mult [n_rows];
	for(i in 1:n_cols) v_mult[,i] = to_array_1d(to_row_vector(v[,i]) .* to_row_vector(r));

	return v_mult;
}

vector pow_vector(vector v, real p){
		vector[rows(v)] v_pow;
		
		for(i in 1:rows(v)) v_pow[i] = v[i] ^ p;
		
		return(v_pow);
	}
	
	real dirichlet_regression_lpdf(vector p, row_vector X, matrix alpha, real phi, real plateau){

		// // Build sum to zero variable
		// int c = cols(alpha);
		// int r = rows(alpha);
		// matrix[r, c]  alpha_ = alpha;
		// alpha_[1,c] = -sum(alpha_[1, 1:(c-1)]);
			//	real  phi_exp= ( 1.0 ./ (phi + 0.0001));

		// Calculate log prob
		return (dirichlet_lpdf(p | softmax( append_row([0]', to_vector(X * alpha))) * exp(phi) + plateau ));
	}

	vector dirichlet_regression_rng( row_vector X, matrix alpha, real phi, real plateau){

		// // Build sum to zero variable
		// int c = cols(alpha);
		// int r = rows(alpha);
		// matrix[r, c]  alpha_ = alpha;
		// alpha_[1,c] = -sum(alpha_[1, 1:(c-1)]);
	//	real  phi_exp= ( 1.0 ./ (phi + 0.0001));

		// Calculate log prob
		return (dirichlet_rng( softmax( append_row([0]', to_vector(X * alpha))) * exp(phi) + plateau ));
	}

real beta_regression_lpdf(vector[] p, matrix X, matrix alpha, real[] phi, real plateau){

		real lp = 0;
		//matrix[num_elements(p[,1]), num_elements(p[1])] mu;
		vector[num_elements(phi)]  phi_exp= ( 1.0 ./ (to_vector(phi )+ 0.0001));

		// Build sum to zero variable
		int c = cols(alpha);
		int r = rows(alpha);
		matrix[r, c+1]  alpha_;
		alpha_[,1:c] = alpha;
		for(rr in 1:r) alpha_[rr,c+1] = -sum(alpha_[rr, 1:c]);
		
	//	print(alpha_);
	//	print(phi);
	//	print(X);
		
		for(j in 1:num_elements(p[,1])) {

			vector[num_elements(p[1])] mu  = softmax( append_row([0]', to_vector(X[j] * alpha)));

	//	print((mu .* phi_exp) +1);


     	lp += beta_lpdf(p[j] | (mu .* phi_exp) +plateau, ((1.0 - mu) .* phi_exp) + plateau);

		}
		return (lp);
	}

vector[] beta_regression_rng( matrix X, matrix alpha, real[] phi, real plateau){

		vector[cols(alpha)+1] p[rows(X)];

		//matrix[num_elements(p[,1]), num_elements(p[1])] mu;
		vector[num_elements(phi)]  phi_exp= ( 1.0 ./ (to_vector(phi )+ 0.0001));

// Build sum to zero variable
		int c = cols(alpha);
		int r = rows(alpha);
		matrix[r, c+1]  alpha_;
		alpha_[,1:c] = alpha;
		for(rr in 1:r) alpha_[rr,c+1] = -sum(alpha_[rr, 1:(c)]);
		
		for(j in 1:num_elements(p[,1])) {


			vector[num_elements(p[1])] mu  = softmax( append_row([0]', to_vector(X[j] * alpha)));

      	 p[j] = to_vector(beta_rng((mu .* phi_exp) +plateau, ((1.0 - mu) .* phi_exp) + plateau));

		}
		return (p);
	}

}
data {
	// shards
	int<lower=1> shards;
	int lv;
	// Reference matrix inference
	int<lower=0> G;
	int<lower=0> GM;
	int<lower=0> S;
	int CL; // counts linear size

	// reference counts
 	int<lower=0> counts_linear[CL] ;

	int G_to_counts_linear[CL] ;
	int S_linear[CL] ;

	// Reference counts per level
	int<lower=0> CL_NA;
	int<lower=0> counts_idx_lv_NA[CL_NA];
	int<lower=0> CL_lv;
	int<lower=0> counts_idx_lv[CL_lv];

	// Priors
	real<upper=0> sigma_slope;
	real<lower=0> sigma_sigma;
	//real sigma_intercept;


	
	// Cell types
 	int<lower=0> Q;
  int<lower=1> n_nodes;
  int<lower=1> ct_in_nodes[n_nodes];
	int<lower=1> n_levels;
  int<lower=1> ct_in_levels[n_levels];

	// reference counts
	int y[Q, GM];
	int max_y;
	matrix[ ct_in_levels[lv], GM] ref;
	
  // Deconvolution
  int<lower=0> G_lv;
  int G_lv_linear[G_lv];

  // Observed counts
  int<lower=0> Y_lv;
	int y_linear_lv[Y_lv];
	int y_linear_S_lv[Y_lv];

	// MPI
	int size_y_linear_MPI[max(shards, 1)];
	int size_y_linear_S_MPI[max(shards, 1)];
	int y_linear_S_MPI[shards,max(size_y_linear_S_MPI)];
	int size_G_linear_MPI[max(shards, 1)];
	int G_linear_MPI[shards,max(size_G_linear_MPI)];
	int y_linear_MPI[shards,max(size_y_linear_MPI)];

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
	vector[G] lambda_log;
  vector[G] sigma_inv_log;

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
	int do_regression;

	// Family
	int fam_dirichlet;

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


  // Proportions
  // lv1
  simplex[ct_in_nodes[1]]  prop_1[Q * (lv == 1)]; // Root

  // lv2
  simplex[ct_in_nodes[2]]  prop_a[Q * (lv == 2)]; // Immune cells childrens

  // lv3
  simplex[ct_in_nodes[3]]  prop_b[Q * (lv == 3)]; // b cells childrens
  simplex[ct_in_nodes[4]]  prop_c[Q * (lv == 3)]; // granulocyte childrens
  simplex[ct_in_nodes[5]]  prop_d[Q * (lv == 3)]; // mono_derived childrens
  simplex[ct_in_nodes[6]]  prop_e[Q * (lv == 3)]; // natural_killer childrens
  simplex[ct_in_nodes[7]]  prop_f[Q * (lv == 3)]; // t_cell childrens

	// lv4
  simplex[ct_in_nodes[8]]  prop_g[Q * (lv == 4)]; // dendritic myeloid childrens
  simplex[ct_in_nodes[9]]  prop_h[Q * (lv == 4)]; // macrophage childrens
  simplex[ct_in_nodes[10]] prop_i[Q * (lv == 4)]; // nk primed
  simplex[ct_in_nodes[11]] prop_l[Q * (lv == 4)]; // CD4 childrens
  simplex[ct_in_nodes[12]] prop_m[Q * (lv == 4)]; // CD8 childrens

	// Dirichlet regression
  // lv1
  matrix[A * (lv == 1) * do_regression,ct_in_nodes[1]-1]  alpha_1; // Root

	// lv2
  matrix[A * (lv == 2) * do_regression,ct_in_nodes[2]-1]  alpha_a; // Immune cells

  // lv3
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[3]-1]  alpha_b; // b cells
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[4]-1]  alpha_c; // granulocyte
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[5]-1]  alpha_d; // mono_derived
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[6]-1]  alpha_e; // natural_killer
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[7]-1]  alpha_f; // t_cell

	// lv4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[8]-1]  alpha_g; // dendritic myeloid
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[9]-1]  alpha_h; // macrophage
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[10]-1]  alpha_i; // NK
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[11]-1] alpha_l; // CD4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[12]-1] alpha_m; // CD8

	real<lower=0> phi[12]; //[fam_dirichlet ? 10 : ct_in_levels[lv]];

	// Unknown population
	row_vector<lower=0, upper = log(max(counts_linear))>[GM] lambda_UFO;
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

 	vector[ct_in_levels[lv]] prop_lv[Q] ;

	// Tree poportion
	vector[ct_in_levels[2]] prop_2[Q * (lv >= 2)];
	vector[ct_in_levels[3]] prop_3[Q * (lv >= 3)];
	vector[ct_in_levels[4]] prop_4[Q * (lv >= 4)];
	
	vector[Q*GM] mu_vector;
	vector[Q*GM] sigma_vector;

	matrix[Q, GM] mu;
	
	real sigma_intercept = 1.5;
	
	// proportion of level 2
	if(lv >= 2)
	prop_2 =
		append_vector_array(
			prop_1_prior[,singles_lv2],
			multiply_by_column( (lv == 2 ? prop_a : prop_a_prior), prop_1_prior[,parents_lv2[1]])
		);

	// proportion of level 3
	if(lv >= 3)
	prop_3 =
		append_vector_array(
			prop_2[,singles_lv3],
			append_vector_array(
				multiply_by_column((lv == 3 ? prop_b : prop_b_prior), prop_2[,parents_lv3[1]]),
				append_vector_array(
					multiply_by_column((lv == 3 ? prop_c : prop_c_prior), prop_2[,parents_lv3[2]]),
					append_vector_array(
						multiply_by_column((lv == 3 ? prop_d : prop_d_prior), prop_2[,parents_lv3[3]]),
						append_vector_array(
							multiply_by_column((lv == 3 ? prop_e : prop_e_prior), prop_2[,parents_lv3[4]]),
							multiply_by_column((lv == 3 ? prop_f : prop_f_prior), prop_2[,parents_lv3[5]])
						)
					)
				)
			)
		);

	// proportion of level 4
	if(lv >= 4)
	prop_4 =
		append_vector_array(
			prop_3[,singles_lv4],
			append_vector_array(
				multiply_by_column(prop_g, prop_3[,parents_lv4[1]]),
				append_vector_array(
					multiply_by_column(prop_h, prop_3[,parents_lv4[2]]),
					append_vector_array(
						multiply_by_column(prop_i, prop_3[,parents_lv4[3]]),
						append_vector_array(
  						multiply_by_column(prop_l, prop_3[,parents_lv4[4]]),
  						multiply_by_column(prop_m, prop_3[,parents_lv4[5]])
  					)
					)
				)
			)
		);
		
	prop_lv	= which(lv, prop_1, prop_2, prop_3, prop_4);

	// Calculate y hat with UFO
	mu =
	
			// Prop matrix
			append_col(	vector_array_to_matrix(prop_lv) * (1-prop_UFO), 	rep_vector(prop_UFO,Q) ) *

			// Expression matrix
			append_row(	ref,	exp(lambda_UFO) );

	// Correct for exposure
	for(q in 1:Q) mu[q] = mu[q] * exposure_rate[q];
	
	// Vectorise 
	mu_vector = to_vector(mu');
	sigma_vector = 1.0 ./ (pow_vector(mu_vector, -0.4) *  exp(sigma_intercept)); //   exp( log(mu_vector)  * -0.4 + sigma_intercept); //;
	
	// Sampling
	to_array_1d(y) ~ neg_binomial_2(mu_vector, sigma_vector);


	// lv 1
  if(lv == 1 && do_regression) {

		//print(X_scaled[,2]);
  	prop_1 ~ beta_regression(X_scaled, alpha_1, phi[1:4], 0.5);
  	 alpha_1[1] ~ normal(0,2);
  	 to_vector( alpha_1[2:] ) ~ student_t(5,  0,  2.5);


  }
	if(lv == 1 && !do_regression) for(q in 1:Q) target += dirichlet_lpdf(prop_1[q] | rep_vector(1, num_elements(prop_1[1])));
	//if(lv > 1)  for(q in 1:Q) target += dirichlet_lpdf(prop_1[q] | prop_1_prior[q]);

	// lv 2
  if(lv == 2 && do_regression) {

  	//prop_a ~ beta_regression(X_scaled, alpha_a, phi[1:6], 1);
  	for(q in 1:Q) prop_a[q] ~ dirichlet_regression( X_scaled[q], alpha_a, phi[1] , 0.05);
  	alpha_a[1] ~ normal(0,2);
  	to_vector( alpha_a[2:] ) ~ student_t(5,  0, 2.5);

  }
	if(lv == 2 && !do_regression) for(q in 1:Q) target += dirichlet_lpdf(prop_a[q] | rep_vector(1, num_elements(prop_a[1])));
	//if(lv > 2)  for(q in 1:Q) target += dirichlet_lpdf(prop_a[q] | prop_a_prior[q]);

	// lv 3
  if(lv == 3 && do_regression){

  		for(q in 1:Q) prop_b[q] ~ dirichlet_regression( X_scaled[q], alpha_b, phi[1] , 1);
  		for(q in 1:Q) prop_c[q] ~ dirichlet_regression( X_scaled[q], alpha_c, phi[2] , 1);
  		for(q in 1:Q) prop_d[q] ~ dirichlet_regression( X_scaled[q], alpha_d, phi[3] , 1);
  		for(q in 1:Q) prop_e[q] ~ dirichlet_regression( X_scaled[q], alpha_e, phi[4] , 1);
  		for(q in 1:Q) prop_f[q] ~ dirichlet_regression( X_scaled[q], alpha_f, phi[5] , 1);

		alpha_b[1] ~  normal(0,2);
  	to_vector( alpha_b[2:] ) ~ student_t(5, 0, 2.5);
		alpha_c[1] ~  normal(0,2);
  	to_vector( alpha_c[2:] ) ~ student_t(5, 0, 2.5);
		alpha_d[1] ~  normal(0,2);
  	to_vector( alpha_d[2:] ) ~ student_t(5, 0, 2.5);
		alpha_e[1] ~  student_t(5, 0, 2.5);
  	to_vector( alpha_e[2:] ) ~ student_t(5, 0, 2.5);
		alpha_f[1] ~  normal(0,2);
  	to_vector( alpha_f[2:] ) ~ student_t(5, 0, 2.5);
  }
  if(lv == 3 && !do_regression) for(q in 1:Q){
  	 target += dirichlet_lpdf(prop_b[q] | rep_vector(1, num_elements(prop_b[1])));
		 target += dirichlet_lpdf(prop_c[q] | rep_vector(1, num_elements(prop_c[1])));
		 target += dirichlet_lpdf(prop_d[q] | rep_vector(1, num_elements(prop_d[1])));
		 target += dirichlet_lpdf(prop_e[q] | rep_vector(1, num_elements(prop_e[1])));
		 target += dirichlet_lpdf(prop_f[q] | rep_vector(1, num_elements(prop_f[1])));
  }

	// lv 4
  if(lv == 4 && do_regression){

  		for(q in 1:Q) prop_g[q] ~ dirichlet_regression( X_scaled[q], alpha_g, phi[1] , 1);
  		for(q in 1:Q) prop_h[q] ~ dirichlet_regression( X_scaled[q], alpha_h, phi[2] , 1);
  		for(q in 1:Q) prop_i[q] ~ dirichlet_regression( X_scaled[q], alpha_i, phi[3] , 1);
  		for(q in 1:Q) prop_l[q] ~ dirichlet_regression( X_scaled[q], alpha_l, phi[4] , 1);
  		for(q in 1:Q) prop_m[q] ~ dirichlet_regression( X_scaled[q], alpha_m, phi[5] , 1);

// else  prop_4 ~ beta_regression(X_scaled, alpha_4, phi);
		alpha_g[1] ~ normal(0,2);
  	to_vector( alpha_g[2:] ) ~ student_t(5, 0, 2.5);
		alpha_h[1] ~ normal(0,2);
  	to_vector( alpha_h[2:] ) ~ student_t(5, 0, 2.5);
		alpha_i[1] ~normal(0,2);
  	to_vector( alpha_i[2:] ) ~ student_t(5, 0, 2.5);
		alpha_l[1] ~ normal(0,2);
  	to_vector( alpha_l[2:] ) ~ student_t(5, 0, 2.5);
		alpha_m[1] ~ normal(0,2);
  	to_vector( alpha_m[2:] ) ~ student_t(5, 0, 2.5);

  }
  if(lv == 4 && !do_regression) for(q in 1:Q){
  	 target += dirichlet_lpdf(prop_g[q] | rep_vector(1, num_elements(prop_g[1])));
		 target += dirichlet_lpdf(prop_h[q] | rep_vector(1, num_elements(prop_h[1])));
		 target += dirichlet_lpdf(prop_i[q] | rep_vector(1, num_elements(prop_i[1])));
		 target += dirichlet_lpdf(prop_l[q] | rep_vector(1, num_elements(prop_l[1])));
		 target += dirichlet_lpdf(prop_m[q] | rep_vector(1, num_elements(prop_m[1])));
  }

	// Dirichlet regression
	phi ~  gamma(1,5);

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
