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
	
	real beta_regression_lpdf(vector[] p, matrix X, matrix beta, vector phi){

		real lp = 0;
    matrix[num_elements(p[1]), num_elements(p[,1])] mu = (X * beta)';
    real buffer;
    
		for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
		for(i in 1:cols(mu)) {
			for(j in 1:rows(mu)) {
			
			buffer = 1.0/fmin(mu[j,i], 1-mu[j,i]) * 1;
			
     	lp += beta_lpdf(p[i,j] | mu[j,i] .* phi[j] * buffer, (1.0 - mu[j,i]) .* phi[j] * buffer );

			// add parashoot
			//lp += beta_lpdf(p[i] | 1.5, 1.5 );
		}}
		return (lp);
	}

	real dirichlet_regression_lpdf(vector[] p, matrix X, matrix beta, real phi, real plateau){
		
		real lp = 0;
    matrix[num_elements(p[1]), num_elements(p[,1])] mu = (X * beta)';
		real buffer;
	
	for(i in 1:cols(mu)){
		mu[,i] = softmax(mu[,i]);
		buffer = 1.0/min(mu[,i]) * plateau;
	
		lp += dirichlet_lpdf(p[i] | mu[,i] * phi * buffer );
	}
	
	return(lp);
	}

  vector Q_sum_to_zero_QR(int N) {
    vector [2*N] Q_r;

    for(i in 1:N) {
      Q_r[i] = -sqrt((N-i)/(N-i+1.0));
      Q_r[i+N] = inv_sqrt((N-i) * (N-i+1));
    }
    return Q_r;
  }

  row_vector sum_to_zero_QR(row_vector x_raw, vector Q_r) {
    int N = num_elements(x_raw) + 1;
    row_vector [N] x;
    real x_aux = 0;

    for(i in 1:N-1){
      x[i] = x_aux + x_raw[i] * Q_r[i];
      x_aux = x_aux + x_raw[i] * Q_r[i+N];
    }
    x[N] = x_aux;
    return x;
  }

	real partial_sum_lpmf(int[] slice_Y,
                        int start, int end,
                        matrix ref,
                        matrix prop,
                        vector exposure_multiplier,
                        real sigma_intercept, 
                        real sigma_slope) {
                          
	matrix[rows(prop), cols(ref)] mu = (prop * ref);
	vector[rows(prop) * cols(ref)] mu_vector;
	for(q in 1:rows(mu)) mu[q] = mu[q] * exposure_multiplier[q];

 mu_vector = to_vector(mu');
	
	return( neg_binomial_2_lupmf(
		slice_Y |
		mu_vector[start:end],
		1.0 ./ (pow_vector(mu_vector[start:end], sigma_slope) *  exp(sigma_intercept))

	));
  

	}
}
data {
	// shards
	int<lower=1> shards;
	int lv;
	
	// Reference matrix inference
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
	vector[spt] prior_survival_time;
	int CIT;
	int columns_idx_including_time[CIT];
	
	// For exposure
	// int nrow_for_exposure;
	// int Q_for_exposure[nrow_for_exposure];
	// int counts_for_exposure[nrow_for_exposure] ;
	// vector[nrow_for_exposure] reference_for_exposure;
	
	// Exposure rate
  vector[Q] exposure_multiplier;

  
}
transformed data{
	real real_data[shards, 0];
	vector[ct_in_nodes[1]*2] Q_r_1 = Q_sum_to_zero_QR(ct_in_nodes[1]);
	vector[ct_in_nodes[2]*2] Q_r_a = Q_sum_to_zero_QR(ct_in_nodes[2]);
	vector[ct_in_nodes[3]*2] Q_r_b = Q_sum_to_zero_QR(ct_in_nodes[3]);
	vector[ct_in_nodes[4]*2] Q_r_c = Q_sum_to_zero_QR(ct_in_nodes[4]);
	vector[ct_in_nodes[5]*2] Q_r_d = Q_sum_to_zero_QR(ct_in_nodes[5]);
	vector[ct_in_nodes[6]*2] Q_r_e = Q_sum_to_zero_QR(ct_in_nodes[6]);
	vector[ct_in_nodes[7]*2] Q_r_f = Q_sum_to_zero_QR(ct_in_nodes[7]);
	vector[ct_in_nodes[8]*2] Q_r_g = Q_sum_to_zero_QR(ct_in_nodes[8]);
	vector[ct_in_nodes[9]*2] Q_r_h = Q_sum_to_zero_QR(ct_in_nodes[9]);
	vector[ct_in_nodes[10]*2] Q_r_i = Q_sum_to_zero_QR(ct_in_nodes[10]);
	vector[ct_in_nodes[11]*2] Q_r_l = Q_sum_to_zero_QR(ct_in_nodes[11]);
	vector[ct_in_nodes[12]*2] Q_r_m = Q_sum_to_zero_QR(ct_in_nodes[12]);
		
		
  real x_raw_sigma = inv_sqrt(1 - inv(ct_in_nodes[lv]));
  
  int grainsize = 1;
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
  matrix[A * (lv == 1) * do_regression,ct_in_nodes[1]-1]  alpha_1_raw; // Root

	// lv2
  matrix[A * (lv == 2) * do_regression,ct_in_nodes[2]-1]  alpha_a_raw; // Immune cells

  // lv3
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[3]-1]  alpha_b_raw; // b cells
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[4]-1]  alpha_c_raw; // granulocyte
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[5]-1]  alpha_d_raw; // mono_derived
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[6]-1]  alpha_e_raw; // natural_killer
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[7]-1]  alpha_f_raw; // t_cell

	// lv4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[8]-1]  alpha_g_raw; // dendritic myeloid
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[9]-1]  alpha_h_raw; // macrophage
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[10]-1] alpha_i_raw; // NK
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[11]-1] alpha_l_raw; // CD4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[12]-1] alpha_m_raw; // CD8

	vector<lower=0>[12] phi; 

	// Unknown population
	row_vector<lower=0, upper = log(max(to_array_1d(y)))>[GM] lambda_UFO;
	vector<lower=0, upper=0.5>[Q] prop_UFO;

	// Censoring
	vector<lower=0>[how_many_cens] unseen;
	real<lower=0> prior_unseen_alpha[how_many_cens > 0];
	real prior_unseen_beta[how_many_cens > 0];
	
// 	 // Local properties of the data
//   vector[Q] exposure_rate;

}
transformed parameters{

  matrix[A * (lv == 1) * do_regression,ct_in_nodes[1]]  alpha_1; // Root
	// lv2
  matrix[A * (lv == 2) * do_regression,ct_in_nodes[2]]  alpha_a; // Immune cells

  // lv3
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[3]]  alpha_b; // b cells
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[4]]  alpha_c; // granulocyte
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[5]]  alpha_d; // mono_derived
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[6]]  alpha_e; // natural_killer
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[7]]  alpha_f; // t_cell

	// lv4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[8]]  alpha_g; // dendritic myeloid
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[9]]  alpha_h; // macrophage
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[10]] alpha_i; // NK
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[11]] alpha_l; // CD4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[12]] alpha_m; // CD8
  
	matrix[Q,A] X_ = X;
	matrix[Q,A] X_scaled = X_;
	
	if(lv == 1) for(a in 1:A)	alpha_1[a] =  sum_to_zero_QR(alpha_1_raw[a], Q_r_1);
	if(lv == 2) for(a in 1:A)	alpha_a[a] =  sum_to_zero_QR(alpha_a_raw[a], Q_r_a);
	if(lv == 3) for(a in 1:A)	alpha_b[a] =  sum_to_zero_QR(alpha_b_raw[a], Q_r_b);
	if(lv == 3) for(a in 1:A)	alpha_c[a] =  sum_to_zero_QR(alpha_c_raw[a], Q_r_c);
	if(lv == 3) for(a in 1:A)	alpha_d[a] =  sum_to_zero_QR(alpha_d_raw[a], Q_r_d);
	if(lv == 3) for(a in 1:A)	alpha_e[a] =  sum_to_zero_QR(alpha_e_raw[a], Q_r_e);
	if(lv == 3) for(a in 1:A)	alpha_f[a] =  sum_to_zero_QR(alpha_f_raw[a], Q_r_f);
	if(lv == 4) for(a in 1:A)	alpha_g[a] =  sum_to_zero_QR(alpha_g_raw[a], Q_r_g);
	if(lv == 4) for(a in 1:A)	alpha_h[a] =  sum_to_zero_QR(alpha_h_raw[a], Q_r_h);
	if(lv == 4) for(a in 1:A)	alpha_i[a] =  sum_to_zero_QR(alpha_i_raw[a], Q_r_i);
	if(lv == 4) for(a in 1:A)	alpha_l[a] =  sum_to_zero_QR(alpha_l_raw[a], Q_r_l);
	if(lv == 4) for(a in 1:A)	alpha_m[a] =  sum_to_zero_QR(alpha_m_raw[a], Q_r_m);

	if(how_many_cens > 0) {
	for(i in 1:CIT)
		for(j in 1:how_many_cens)
			if(X_[which_cens[j],columns_idx_including_time[i]] > 0)
				X_[which_cens[j],columns_idx_including_time[i]] = sqrt(X_[which_cens[j],columns_idx_including_time[i]]^2 + unseen[j]^2);

	
	
	for(i in 1:CIT)	
			 X_scaled[,columns_idx_including_time[i]] = (X_scaled[,columns_idx_including_time[i]] - mean(X_scaled[,columns_idx_including_time[i]])) / sd(X_scaled[,columns_idx_including_time[i]]);


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
	
	real sigma_intercept = 1.3420415;

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


 target += reduce_sum(
  	partial_sum_lupmf,  to_array_1d(y), grainsize,
  	append_row(	ref,	exp(lambda_UFO) ),
  	append_col(	multiply_matrix_by_column( prop_lv,  (1-prop_UFO) ), 	prop_UFO ),
  	exposure_multiplier,
  	sigma_intercept,
  	-0.4
  );


	// lv 1
  if(lv == 1 && do_regression) {

  	 prop_1 ~ beta_regression(X_scaled, alpha_1, exp(phi[1:4]));
  	 alpha_1_raw[1] ~ normal(0, 2);
  	 if(A > 1) to_vector( alpha_1_raw[2:] ) ~ normal(0, 2);


  }
	if(lv == 1 && !do_regression) for(q in 1:Q) target += dirichlet_lpdf(prop_1[q] | rep_vector(1, num_elements(prop_1[1])));

	// lv 2
  if(lv == 2 && do_regression) {

  	//prop_a ~ beta_regression(X_scaled, alpha_a, phi[1:6], 1);
  	prop_a ~ dirichlet_regression( X_scaled, alpha_a, exp(phi[1]) , 0.5);
  	alpha_a_raw[1] ~ normal(0,2);
  	if(A > 1)  to_vector( alpha_a_raw[2:] ) ~ normal(0, 2);

  }
	if(lv == 2 && !do_regression) for(q in 1:Q) target += dirichlet_lpdf(prop_a[q] | rep_vector(1, num_elements(prop_a[1])));

	// lv 3
  if(lv == 3 && do_regression){

  		prop_b ~ dirichlet_regression( X_scaled, alpha_b,  exp(phi[1]) , 1);
  		prop_c ~ dirichlet_regression( X_scaled, alpha_c, exp(phi[2]) , 1);
  		prop_d ~ dirichlet_regression( X_scaled, alpha_d,  exp(phi[3]) , 1);
  		prop_e ~ dirichlet_regression( X_scaled, alpha_e,  exp(phi[4]), 1);
  		prop_f ~ dirichlet_regression( X_scaled, alpha_f,  exp(phi[5]), 1);

		to_vector(alpha_b_raw) ~  normal(0,1);
		to_vector(alpha_c_raw) ~  normal(0,1);
		to_vector(alpha_d_raw) ~  normal(0,1);
		to_vector(alpha_e_raw) ~  normal(0,1);
		to_vector(alpha_f_raw) ~  normal(0,1);
					
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

  		 prop_g ~ dirichlet_regression( X_scaled, alpha_g,  exp(phi[1]) , 1);
  		 prop_h ~ dirichlet_regression( X_scaled, alpha_h,  exp(phi[2]) , 1);
  		 prop_i ~ dirichlet_regression( X_scaled, alpha_i,  exp(phi[3]) , 1);
  		 prop_l ~ dirichlet_regression( X_scaled, alpha_l,  exp(phi[4]), 1);
  		 prop_m ~ dirichlet_regression( X_scaled, alpha_m,  exp(phi[5]) , 1);

// else  prop_4 ~ beta_regression(X_scaled, alpha_4, phi);
		to_vector(alpha_g_raw) ~  normal(0,1);
		to_vector(alpha_h_raw) ~  normal(0,1);
		to_vector(alpha_i_raw) ~  normal(0,1);
		to_vector(alpha_l_raw) ~  normal(0,1);
		to_vector(alpha_m_raw) ~  normal(0,1);

  }
  if(lv == 4 && !do_regression) for(q in 1:Q){
  	 target += dirichlet_lpdf(prop_g[q] | rep_vector(1, num_elements(prop_g[1])));
		 target += dirichlet_lpdf(prop_h[q] | rep_vector(1, num_elements(prop_h[1])));
		 target += dirichlet_lpdf(prop_i[q] | rep_vector(1, num_elements(prop_i[1])));
		 target += dirichlet_lpdf(prop_l[q] | rep_vector(1, num_elements(prop_l[1])));
		 target += dirichlet_lpdf(prop_m[q] | rep_vector(1, num_elements(prop_m[1])));
  }

	phi ~ normal(5,2);

	// lambda UFO
	for(i in 1:shards) lambda_UFO[i] ~ skew_normal(6.2, 3.3, -2.7);
	target += beta_lpdf(prop_UFO | 1.001, 20);

	// Censoring

	if(how_many_cens > 0){
		

		// unseen
		// unseen ~ gamma(1,2);
		X_[which_cens,2] ~ gamma( prior_unseen_alpha[1], prior_unseen_beta[1]);

		// Priors
		//target += gamma_lpdf(X[which_not_cens,2] | prior_unseen_alpha[1], prior_unseen_beta[1]);
	 	target += gamma_lccdf(	X_[which_cens,2] | prior_unseen_alpha[1], prior_unseen_beta[1]);
	 	
	 	// Hyperprior
	 	prior_survival_time ~ gamma( prior_unseen_alpha[1], prior_unseen_beta[1]);
	 	
	 	prior_unseen_beta[1] ~ student_t(3, 6, 10);
  	prior_unseen_alpha[1] ~  gamma(0.01, 0.01);

	}
}
