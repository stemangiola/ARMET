functions{

	
matrix vector_array_to_matrix(vector[] x) {
		matrix[size(x), rows(x[1])] y;
		for (m in 1:size(x))
		  y[m] = x[m]';
		return y;
}

matrix which(int x, vector[] a, matrix b, matrix c, matrix d){
		if(x == 1) return(vector_array_to_matrix(a));
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

matrix multiply_by_column(vector[] v, vector r){
	int n_rows = num_elements(v[,1]);
	int n_cols = num_elements(v[1]);

	matrix[n_rows,n_cols] v_mult ;
	for(i in 1:n_rows) v_mult[i] = to_row_vector(v[i] * r[i]);

	return v_mult;
}

matrix multiply_matrix_by_column(matrix v, vector r){
	int n_rows = rows(v);
	int n_cols = cols(v);

	matrix[n_rows,n_cols] v_mult ;
	for(i in 1:n_rows) v_mult[i] = v[i] * r[i];

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
	
	// SUM TO ZERO FRAMEWORK
	vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }

  vector sum_to_zero_QR(vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }
  
  vector[] rate_to_prop(vector[] v, vector Q_r){
  		int R = num_elements(v[,1]);
			int C = num_elements(v[1]);
	
  	vector[C+1] prop[R];
  	
  	for(r in 1:R)
  		prop[r] = softmax(sum_to_zero_QR(v[r], Q_r));
  		
  	return(prop);
  }
  

}
data {
	// shards
	int<lower=1> shards;
	int lv;
	
	// Reference matrix inference
	int<lower=0> G;
	int<lower=0> GM;

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
	int ct_in_ancestor_level;
	matrix[ ct_in_levels[lv], GM] ref;
	matrix[ Q, ct_in_ancestor_level] prior_prop;
	
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

	// Censoring
	int how_many_cens;
	int which_cens[how_many_cens];
	int which_not_cens[Q-how_many_cens];
	real<lower=0> max_unseen;
	int spt;
	real prior_survival_time[spt];

  // Local properties of the data
  vector[Q] exposure_rate;
  
}
transformed data{

	vector[2*ct_in_nodes[1]] Q_r_1 = Q_sum_to_zero_QR(ct_in_nodes[1]);
	vector[2*ct_in_nodes[2]] Q_r_a = Q_sum_to_zero_QR(ct_in_nodes[2]);
	vector[2*ct_in_nodes[3]] Q_r_b = Q_sum_to_zero_QR(ct_in_nodes[3]);
	vector[2*ct_in_nodes[4]] Q_r_c = Q_sum_to_zero_QR(ct_in_nodes[4]);
	vector[2*ct_in_nodes[5]] Q_r_d = Q_sum_to_zero_QR(ct_in_nodes[5]);
	vector[2*ct_in_nodes[6]] Q_r_e = Q_sum_to_zero_QR(ct_in_nodes[6]);
	vector[2*ct_in_nodes[7]] Q_r_f = Q_sum_to_zero_QR(ct_in_nodes[7]);
	vector[2*ct_in_nodes[8]] Q_r_g = Q_sum_to_zero_QR(ct_in_nodes[8]);
	vector[2*ct_in_nodes[9]] Q_r_h = Q_sum_to_zero_QR(ct_in_nodes[9]);
	vector[2*ct_in_nodes[10]] Q_r_i = Q_sum_to_zero_QR(ct_in_nodes[10]);
	vector[2*ct_in_nodes[11]] Q_r_l = Q_sum_to_zero_QR(ct_in_nodes[11]);
	vector[2*ct_in_nodes[12]] Q_r_m = Q_sum_to_zero_QR(ct_in_nodes[12]);

}
parameters {


  // Proportions rates
  // lv1
  simplex[ct_in_nodes[1]-1]  rate_1[Q * (lv == 1)]; // Root

  // lv2
  simplex[ct_in_nodes[2]-1]  rate_a[Q * (lv == 2)]; // Immune cells childrens

  // lv3
  simplex[ct_in_nodes[3]-1]  rate_b[Q * (lv == 3)]; // b cells childrens
  simplex[ct_in_nodes[4]-1]  rate_c[Q * (lv == 3)]; // granulocyte childrens
  simplex[ct_in_nodes[5]-1]  rate_d[Q * (lv == 3)]; // mono_derived childrens
  simplex[ct_in_nodes[6]-1]  rate_e[Q * (lv == 3)]; // natural_killer childrens
  simplex[ct_in_nodes[7]-1]  rate_f[Q * (lv == 3)]; // t_cell childrens

	// lv4
  simplex[ct_in_nodes[8]-1]  rate_g[Q * (lv == 4)]; // dendritic myeloid childrens
  simplex[ct_in_nodes[9]-1]  rate_h[Q * (lv == 4)]; // macrophage childrens
  simplex[ct_in_nodes[10]-1] rate_i[Q * (lv == 4)]; // nk primed
  simplex[ct_in_nodes[11]-1] rate_l[Q * (lv == 4)]; // CD4 childrens
  simplex[ct_in_nodes[12]-1] rate_m[Q * (lv == 4)]; // CD8 childrens


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
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[10]-1] alpha_i; // NK
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[11]-1] alpha_l; // CD4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[12]-1] alpha_m; // CD8

	real<lower=0> phi[12]; 

	// Unknown population
	row_vector<lower=0, upper = log(max(to_array_1d(y)))>[GM] lambda_UFO;
	vector<lower=0, upper=1>[Q] prop_UFO;

	// Censoring
	vector<lower=0>[how_many_cens] unseen;
	real<lower=0> prior_unseen_alpha[how_many_cens > 0];
	real prior_unseen_beta[how_many_cens > 0];

}
transformed parameters{

	// Proportions
  // lv1
  
  vector[ct_in_nodes[1]]  prop_1[Q * (lv == 1)] = rate_to_prop(rate_1, Q_r_1); // Root

  // lv2
  vector[ct_in_nodes[2]]  prop_a[Q * (lv == 2)] = rate_to_prop(rate_a, Q_r_a); // Immune cells childrens

  // lv3
  vector[ct_in_nodes[3]]  prop_b[Q * (lv == 3)] = rate_to_prop(rate_b, Q_r_b); // b cells childrens
  vector[ct_in_nodes[4]]  prop_c[Q * (lv == 3)] = rate_to_prop(rate_c, Q_r_c); // granulocyte childrens
  vector[ct_in_nodes[5]]  prop_d[Q * (lv == 3)] = rate_to_prop(rate_d, Q_r_d); // mono_derived childrens
  vector[ct_in_nodes[6]]  prop_e[Q * (lv == 3)] = rate_to_prop(rate_e, Q_r_e); // natural_killer childrens
  vector[ct_in_nodes[7]]  prop_f[Q * (lv == 3)] = rate_to_prop(rate_f, Q_r_f); // t_cell childrens

	// lv4
  vector[ct_in_nodes[8]]  prop_g[Q * (lv == 4)] = rate_to_prop(rate_g, Q_r_g); // dendritic myeloid childrens
  vector[ct_in_nodes[9]]  prop_h[Q * (lv == 4)] = rate_to_prop(rate_h, Q_r_h); // macrophage childrens
  vector[ct_in_nodes[10]] prop_i[Q * (lv == 4)] = rate_to_prop(rate_i, Q_r_i); // nk primed
  vector[ct_in_nodes[11]] prop_l[Q * (lv == 4)] = rate_to_prop(rate_l, Q_r_l); // CD4 childrens
  vector[ct_in_nodes[12]] prop_m[Q * (lv == 4)] = rate_to_prop(rate_m, Q_r_m); // CD8 childrens
  
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

 	matrix[Q, ct_in_levels[lv]] prop_lv ;

	// Tree poportion
	matrix[Q * (lv >= 2), ct_in_levels[2]] prop_2;
	matrix[Q * (lv >= 2), ct_in_levels[3]] prop_3;
	matrix[Q * (lv >= 2), ct_in_levels[4]] prop_4;
	
	vector[Q*GM] mu_vector;
	vector[Q*GM] sigma_vector;

	matrix[Q, GM] mu;
	
	real sigma_intercept = 1.5;
	
  // Compensation QR sum-to-zero
  if(lv ==1) for(q in 1:Q) rate_1[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[1])));
  	
  if(lv ==2) for(q in 1:Q) rate_a[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[2])));

  if(lv ==3) {
  	for(q in 1:Q) rate_b[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[3])));
  	for(q in 1:Q) rate_c[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[4])));
  	for(q in 1:Q) rate_d[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[5])));
  	for(q in 1:Q) rate_e[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[6])));
  	for(q in 1:Q) rate_f[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[7])));
  }

  if(lv ==4) {
  	for(q in 1:Q) rate_g[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[8])));
  	for(q in 1:Q) rate_h[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[9])));
  	for(q in 1:Q) rate_i[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[10])));
  	for(q in 1:Q) rate_l[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[11])));
  	for(q in 1:Q) rate_m[q] ~ normal(0, inv_sqrt(1 - inv(ct_in_nodes[12])));
  }
  
	// proportion of level 2
	if(lv == 2)
	prop_2 =
		append_col(
			prior_prop[,singles_lv2],
			multiply_by_column( prop_a , prior_prop[,parents_lv2[1]])
		);

	// proportion of level 3
	if(lv == 3)
	prop_3 =
		append_col(
			prior_prop[,singles_lv3],
			append_col(
				multiply_by_column(prop_b , prior_prop[,parents_lv3[1]]),
				append_col(
					multiply_by_column(prop_c , prior_prop[,parents_lv3[2]]),
					append_col(
						multiply_by_column(prop_d , prior_prop[,parents_lv3[3]]),
						append_col(
							multiply_by_column( prop_e, prior_prop[,parents_lv3[4]]),
							multiply_by_column(prop_f , prior_prop[,parents_lv3[5]])
						)
					)
				)
			)
		);

	// proportion of level 4
	if(lv == 4)
	prop_4 =
		append_col(
			prior_prop[,singles_lv4],
			append_col(
				multiply_by_column(prop_g, prior_prop[,parents_lv4[1]]),
				append_col(
					multiply_by_column(prop_h, prior_prop[,parents_lv4[2]]),
					append_col(
						multiply_by_column(prop_i, prior_prop[,parents_lv4[3]]),
						append_col(
  						multiply_by_column(prop_l, prior_prop[,parents_lv4[4]]),
  						multiply_by_column(prop_m, prior_prop[,parents_lv4[5]])
  					)
					)
				)
			)
		);
		
	prop_lv	= which(lv, prop_1, prop_2, prop_3, prop_4);

	// Calculate y hat with UFO
	mu =
	
			// Prop matrix
			append_col(	multiply_matrix_by_column( prop_lv,  (1-prop_UFO) ), 	prop_UFO ) *

			// Expression matrix
			append_row(	ref,	exp(lambda_UFO) );

	// Correct for exposure
	for(q in 1:Q) mu[q] = mu[q] * 1.0/exposure_rate[q];
	
	// Sampling
		
	// Vectorise 
	mu_vector = to_vector(mu');
	sigma_vector = 1.0 ./ (pow_vector(mu_vector, -0.4) *  exp(sigma_intercept)); //   exp( log(mu_vector)  * -0.4 + sigma_intercept); //;
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
	target += beta_lpdf(prop_UFO | 1.001, 20);

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
