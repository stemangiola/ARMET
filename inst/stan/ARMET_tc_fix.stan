functions{

	matrix vector_array_to_matrix(vector[] x) {
			matrix[size(x), rows(x[1])] y;
			for (m in 1:size(x))
			  y[m] = x[m]';
			return y; 
	}
	
	matrix vector_array_to_matrix_transpose(vector[] x) {
			matrix[rows(x[1]), size(x)] y;
			for (m in 1:cols(y))
			  y[,m] = x[m];
			return y; 
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
	
	matrix multiply_matrix_by_row(matrix v, row_vector r){
		int n_rows = rows(v);
		int n_cols = cols(v);
	
		matrix[n_rows,n_cols] v_mult ;
		for(i in 1:n_cols) v_mult[,i] = v[,i] * r[i];
	
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
    vector[rows(mu)] buffer; // one buffer for each phi, j is a cell type
    real sample_with_min_proportion;
    
    // Softmax
		for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
		
		// Buffer, to make the alpha and beta > 1, and so unimodal
		for(j in 1:rows(mu)) {
			sample_with_min_proportion = min(mu[j,]);
			buffer[j] = 1.0/fmin(sample_with_min_proportion, 1-sample_with_min_proportion) * 1;
		}
		
		for(i in 1:cols(mu)) {
			for(j in 1:rows(mu)) {

     	lp += beta_lpdf(p[i,j] | mu[j,i] .* phi[j] , (1.0 - mu[j,i]) .* phi[j] );

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
  
  matrix get_A_qr(int K){
  		matrix[K, K] A = diag_matrix(rep_vector(1,K));
  
		 for (i in 1:K-1) A[K,i] = -1;
  		A[K,K] = 0;
  	
  	return qr_Q(A)[ , 1:(K-1)];
  
  }
  
   vector sum_to_zero_QR_2(matrix A_qr, vector beta_raw) {
    
    return A_qr * beta_raw;
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
	
		real partial_sum2_lpmf(int[] slice_Y,
                        int start, int end,
                        matrix ref,
                        matrix prop,
                        vector exposure_multiplier,
                        real sigma_intercept, 
                        real sigma_slope) {
                         
	matrix[rows(ref), cols(prop)] mu = (ref * prop );
	vector[cols(prop) * rows(ref)] mu_vector;
	for(q in 1:cols(mu)) mu[,q] = mu[,q] * exposure_multiplier[q];

 mu_vector = to_vector(mu)[start:end];
	
	return( neg_binomial_2_lupmf(
		slice_Y |
		mu_vector,
		1.0 ./ (pow_vector(mu_vector, sigma_slope) *  exp(sigma_intercept))

	));
  

	}
	
}
data {
	// shards
	int<lower=1> shards;

		
	// Reference matrix inference
	int<lower=0> GM;

	// Priors
	real<upper=0> sigma_slope;
	real<lower=0> sigma_sigma;
	//real sigma_intercept;
	
	// Cell types
 	int<lower=0> Q;
 	int<lower=2> number_of_cell_types;
 	
	// reference counts
	int y[Q, GM];
	int max_y;
	matrix[ number_of_cell_types, GM] ref;

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
  
  int<lower=0, upper=1> use_data;

  
}
transformed data{
	int lv = 1;
	matrix[GM, number_of_cell_types]  ref_t = ref';
	int y_array[Q * GM] = to_array_1d(y);
	vector[number_of_cell_types*2] Q_r_1 = Q_sum_to_zero_QR(number_of_cell_types);
	matrix[number_of_cell_types, number_of_cell_types-1] A_qr = get_A_qr(number_of_cell_types);

  real x_raw_sigma = inv_sqrt(1 - inv(number_of_cell_types));
  
  int grainsize = 1;
}
parameters {


  // Proportions
  // lv1
  matrix[number_of_cell_types-1, Q]  prop_1_raw; // Root

	// Dirichlet regression
  // lv1
  matrix[A  * do_regression,number_of_cell_types-1]  alpha_1_raw; // Root

	vector<lower=1>[number_of_cell_types] phi; 

	// Unknown population
	vector<lower=0, upper = log(max(to_array_1d(y)))>[GM] lambda_UFO;
	row_vector<lower=0, upper=0.5>[Q] prop_UFO;

	// Censoring
	vector<lower=0>[how_many_cens] unseen;
	real<lower=0> prior_unseen_alpha[how_many_cens > 0];
	real prior_unseen_beta[how_many_cens > 0];
	
// 	 // Local properties of the data
//   vector[Q] exposure_rate;

}
transformed parameters{

  matrix[A  * do_regression,number_of_cell_types]  alpha_1; // Root

	matrix[Q,A] X_ = X;
	matrix[Q,A] X_scaled = X;
	matrix[number_of_cell_types, Q]  prop_1;
	for(q in 1:Q) prop_1[,q] =  softmax( sum_to_zero_QR_2( A_qr, prop_1_raw[,q] )) ;
	
	
	for(a in 1:A)	alpha_1[a] =  sum_to_zero_QR(alpha_1_raw[a], Q_r_1);
	

	if(how_many_cens > 0) {
	for(i in 1:CIT)
		for(j in 1:how_many_cens)
			if(X_[which_cens[j],columns_idx_including_time[i]] > 0)
				X_[which_cens[j],columns_idx_including_time[i]] = sqrt(X_[which_cens[j],columns_idx_including_time[i]]^2 + unseen[j]^2);

	
	X_scaled = X_;

	for(i in 1:CIT)	
			 X_scaled[,columns_idx_including_time[i]] = (X_scaled[,columns_idx_including_time[i]] - mean(X_scaled[,columns_idx_including_time[i]])) / sd(X_scaled[,columns_idx_including_time[i]]);


}

}
model {

	real sigma_intercept = 1.3420415;
		
	if(use_data==1)
			target += reduce_sum(
		  	partial_sum2_lupmf,  y_array, grainsize,
		  	append_col(	ref_t,	exp(lambda_UFO) ),
		  	append_row(	multiply_matrix_by_row( prop_1,  (1-prop_UFO) ), 	prop_UFO ),
		  	exposure_multiplier,
		  	sigma_intercept,
		  	-0.4
		  );


	// lv 1
  if(do_regression) {

  	 for(c in 1:(number_of_cell_types-1)) prop_1_raw[c] ~  normal(X_scaled * alpha_1_raw[,c] * 0.5, phi[c]);
  	 to_vector(alpha_1_raw) ~ normal(0, x_raw_sigma);
  	 //if(A > 1) to_vector( alpha_1_raw[2:] ) ~ normal(0, 0.5);


  }
// 	else for(q in 1:Q) prop_1[q] ~ dirichlet(rep_vector(1, num_elements(prop_1[1])));


	phi ~ normal(0,1);

	// lambda UFO
	lambda_UFO ~ skew_normal(6.2, 3.3, -2.7);
	prop_UFO ~ beta( 1.001, 20);

	// Censoring

	if(how_many_cens > 0){
		
		// Priors unseen
		target += gamma_lpdf(X_[which_not_cens,2] | prior_unseen_alpha[1], prior_unseen_beta[1]);
	 	target += gamma_lccdf(X_[which_cens,2] | prior_unseen_alpha[1], prior_unseen_beta[1]);
	 	
	 	// Hyperprior
	 	prior_survival_time ~ gamma( prior_unseen_alpha[1], prior_unseen_beta[1]);
	 	
	 	prior_unseen_beta[1] ~ student_t(3, 6, 10);
  	prior_unseen_alpha[1] ~  gamma(0.01, 0.01);

	}
}
