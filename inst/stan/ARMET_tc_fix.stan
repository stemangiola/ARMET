functions{

  vector[] get_mean_prop(matrix X, matrix alpha){

	  	vector[cols(alpha)] mu[rows(X)];

		  for(j in 1:num_elements(mu[,1])) 
	  		mu[j]  = softmax( to_vector(X[j] * alpha));
  
  	return(mu);
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
  
  matrix get_A_qr(int K){
  		matrix[K, K] A = diag_matrix(rep_vector(1,K));
  
		 for (i in 1:K-1) A[K,i] = -1;
  		A[K,K] = 0;
  	
  	return qr_Q(A)[ , 1:(K-1)];
  
  }
  
   vector sum_to_zero_QR_2(matrix A_qr, vector beta_raw) {
    
    return A_qr * beta_raw;
  }


	
		real partial_sum2_lpmf(int[] slice_Y,
                        int start, int end,
                        matrix ref,
                        matrix prop,
                        vector exposure_multiplier,
                        real sigma_intercept, 
                        real sigma_slope) {
                         
	matrix[rows(ref), cols(prop)] mu = (ref * prop );
	vector[end-start+1] mu_vector;
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
	//real<upper=0> sigma_slope;
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
  
  real lambda_UFO_mean = mean(log1p(to_array_1d(ref)));
  real lambda_UFO_sd = sd(log1p(to_array_1d(ref)));
  
}
parameters {


  // Proportions
  matrix[number_of_cell_types-1, Q]  prop_1_raw_raw; // Root

	// Dirichlet regression
  matrix[A  * do_regression,number_of_cell_types-1]  alpha_1_raw; // Root

	vector<lower=0>[number_of_cell_types] phi; 

	real sigma_intercept;
	real<upper=0> sigma_slope;
	
	// Unknown population
	vector<lower=0, upper = log(max(to_array_1d(y)))>[GM] lambda_UFO;
	row_vector<lower=0, upper=0.5>[Q] prop_UFO;

	// Censoring
	vector<lower=0>[how_many_cens] unseen;
	real<lower=0> prior_unseen_alpha[how_many_cens > 0];
	real prior_unseen_beta[how_many_cens > 0];



}
transformed parameters{


	matrix[Q,A] X_ = X;
	matrix[Q,A] X_scaled = X;
	matrix[number_of_cell_types, Q]  prop_1;
	matrix[number_of_cell_types-1, Q]  prop_1_raw; // Root

	for(c in 1:(number_of_cell_types-1)) prop_1_raw[c] = to_row_vector(X_scaled * alpha_1_raw[,c]) + phi[c] * prop_1_raw_raw[c];
	
	for(q in 1:Q) prop_1[,q] =  softmax( sum_to_zero_QR(  prop_1_raw[,q], Q_r_1 )) ;
	

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


		
	if(use_data==1){
		
			matrix[GM, Q] mu = // matrix G x Q
				append_col(	ref_t,	exp(lambda_UFO) ) * // MU
				append_row(	multiply_matrix_by_row( prop_1,  (1-prop_UFO) ), 	prop_UFO ); // PROP
				
			for(q in 1:Q) mu[,q] = mu[,q] * exposure_multiplier[q];
			
			vector[GM * Q] mu_vector = to_vector(mu);

			y_array ~  neg_binomial_2(mu_vector,	1.0 ./ (pow_vector(mu_vector, sigma_slope) *  exp(sigma_intercept)));
	}

	// Prior
  to_vector(prop_1_raw_raw) ~  std_normal();

	// Hyper Prior
	to_vector(alpha_1_raw) ~ normal(0,3);
	phi ~ gamma(1.1, 5);
	sigma_intercept ~ normal(0,1);
	sigma_slope ~ normal(0,0.5);
	

	// lambda UFO
	lambda_UFO ~ normal( lambda_UFO_mean , lambda_UFO_sd);
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

generated quantities{

	matrix[A  * do_regression,number_of_cell_types]  alpha_1; // Root
  matrix[number_of_cell_types-1, Q]  prop_1_raw_rng; 
	matrix[number_of_cell_types, Q]  prop_1_rng;
	vector[number_of_cell_types]  mu_1_rng[Q  ]; 
	
	for(c in 1:(number_of_cell_types-1)) prop_1_raw_rng[c] =  to_row_vector(normal_rng(X_scaled * alpha_1_raw[,c] * 0.5, phi[c]));
	for(q in 1:Q) prop_1_rng[,q] =  softmax( sum_to_zero_QR( prop_1_raw_rng[,q], Q_r_1 )) ;
	for(a in 1:A)	alpha_1[a] =  to_row_vector(sum_to_zero_QR(  to_vector(alpha_1_raw[a]), Q_r_1 ));
	
	mu_1_rng = get_mean_prop(X_scaled, alpha_1);

 
}
