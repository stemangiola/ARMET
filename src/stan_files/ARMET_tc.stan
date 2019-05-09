functions{

		matrix vector_array_to_transposed_matrix(vector[] x) {
		  matrix[rows(x[1]), size(x)] y;
		  for (m in 1:size(x))
		    y[,m] = x[m];
		  return y;
		}

		vector vector_array_to_vector(vector[] x) {
			// This operation is column major
		  vector[rows(x[1]) * size(x)] y;
		  for (m in 1:rows(x[1]))
		    y[(m-1)*size(x)+1:(m*size(x))] = to_vector(x[,m]);
		  return y;
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


		vector sum_NB(vector lambda, vector sigma, int[] matrix_dim, vector[] prop){

			int Q = size(prop);
			matrix[matrix_dim[1], matrix_dim[2]] lambda_mat = to_matrix(lambda, matrix_dim[1], matrix_dim[2]); // Bug prone
			matrix[rows(prop[1]),Q] prop_mat = vector_array_to_transposed_matrix(prop);
			matrix[matrix_dim[1] , Q] lambda_sum = lambda_mat * prop_mat;
			matrix[matrix_dim[1], Q] sigma_sum =
				square(lambda_sum) ./
				((
					square(lambda_mat) ./
					to_matrix(sigma, matrix_dim[1], matrix_dim[2]) //sigma_mat
				) * square(prop_mat)) ;

				return(append_row( to_vector(lambda_sum), to_vector(sigma_sum)));
		}

  	real exp_gamma_meanSd_lpdf(vector x_log, real m_log, real s){

      // This function is the  probability of the log gamma function
      // in case you have data that is aleady in log form

      real m = exp(m_log);
      real v = m + square(m) * s;
      real a = square(m) / v;
      real b = m / v;

  		vector[rows(x_log)] jacob = x_log; //jacobian
  		real norm_constant = a * log(b) -lgamma(a);
  		real a_minus_1 = a-1;
  		return sum( jacob ) + norm_constant * rows(x_log) + sum(  x_log * a_minus_1 - exp(x_log) * b ) ;

  	}

  	real exp_gamma_meanSd_rng(real m_log, real s){
  	  // This function takes care of the two prior choice
  	  // without complicating too much the model itself

      real m = exp(m_log);
      real v = m + square(m) * s;
      real a = square(m) / v;
      real b = m / v;

  	  return log(gamma_rng(a, b));
  	}

	vector lp_reduce( vector global_parameters , vector local_parameters , real[] xr , int[] xi ) {
	 	int M = xi[1];
	 	int N = xi[2];
	 	int S = xi[3];
	 	int G_per_shard = xi[4];
	 	int symbol_end[M+1] = xi[(4+1):(4+1+M)];
	 	int sample_idx[N] = xi[(4+1+M+1):(4+1+M+1+N-1)];
	 	int counts[N] = xi[(4+1+M+1+N):size(xi)];


	 	vector[G_per_shard] lambda_MPI = local_parameters[1:G_per_shard];
	 	vector[G_per_shard] sigma_MPI = local_parameters[(G_per_shard+1):(G_per_shard*2)];
	 	vector[S] exposure_rate = local_parameters[((M*2)+1):rows(local_parameters)];

	 	vector[G_per_shard] lp;


	 for(g in 1:G_per_shard){
	 	lp[g] =  neg_binomial_2_log_lpmf(
	 	  counts[(symbol_end[g]+1):symbol_end[g+1]] |
	 	  exposure_rate[sample_idx[(symbol_end[g]+1):symbol_end[g+1]]] +
	 	  lambda_MPI[g],
	 	  sigma_MPI[g]
	 	 );
	 }


    return [sum(lp)]';
  }

		vector sum_NB_MPI(matrix lambda_mat, matrix sigma_mat, matrix prop_mat){

			matrix[rows(lambda_mat) , cols(prop_mat)] lambda_sum = lambda_mat * prop_mat;
			matrix[rows(sigma_mat), cols(prop_mat)] sigma_sum =
				square(lambda_sum) ./
				((
					square(lambda_mat) ./
					sigma_mat
				) * square(prop_mat)) ;

				return(append_row( to_vector(lambda_sum), to_vector(sigma_sum)));
		}

	vector sum_reduce( vector global_parameters , vector local_parameters , real[] xr , int[] xi ) {

		int ct_in_levels = xi[1];
		int Q = xi[2];
		int S = xi[3];
		int y_MPI_symbol_per_shard = xi[4];
		int y_MPI_G_per_shard = xi[5];
		int y_MPI_N_per_shard = xi[6];
	 	int counts[y_MPI_N_per_shard] = xi[6+1: 6+y_MPI_N_per_shard];

	 	vector[y_MPI_G_per_shard] lambda_MPI = local_parameters[1:y_MPI_G_per_shard];
	 	vector[y_MPI_G_per_shard] sigma_MPI = local_parameters[(y_MPI_G_per_shard+1):(y_MPI_G_per_shard*2)];
	 	vector[Q] exposure_rate = local_parameters[((y_MPI_G_per_shard*2)+1):((y_MPI_G_per_shard*2)) + Q];
	 	vector[Q * ct_in_levels] prop = local_parameters[((y_MPI_G_per_shard*2)) + Q+1	:	(((y_MPI_G_per_shard*2)) + Q) + (Q * ct_in_levels)];


		vector[y_MPI_N_per_shard * 2] my_sum = sum_NB_MPI(
			to_matrix( lambda_MPI, y_MPI_symbol_per_shard, ct_in_levels, 0), // Row major
			to_matrix( sigma_MPI,  y_MPI_symbol_per_shard, ct_in_levels, 0), // Row major
			to_matrix( prop, Q, ct_in_levels)'
		);


// print( to_matrix( lambda_MPI, y_MPI_symbol_per_shard, ct_in_levels, 0));
// print(to_matrix( sigma_MPI,  y_MPI_symbol_per_shard, ct_in_levels, 0));
// print(to_matrix( prop, Q, ct_in_levels)');
//
// print(my_sum[1:y_MPI_N_per_shard]);
// print(my_sum[(y_MPI_N_per_shard+1):(y_MPI_N_per_shard*2)]);
// print(to_vector(rep_matrix(to_row_vector(exp(exposure_rate)), y_MPI_symbol_per_shard)));

		// Vecotrised sampling
		return [
			neg_binomial_2_lpmf(
				counts |
				my_sum[1:y_MPI_N_per_shard] .*  to_vector(rep_matrix(to_row_vector(exp(exposure_rate)), y_MPI_symbol_per_shard)),
				my_sum[(y_MPI_N_per_shard+1):(y_MPI_N_per_shard*2)]
			)]';
  }
}
data {

	// Reference matrix inference
  int<lower=0> N;
  int<lower=0> M;
	int<lower=0> G;
	int<lower=0> S;
  int n_shards;
	int<lower=0> counts[n_shards, N];
	int<lower=0> symbol_end[n_shards, M+1];
	int<lower=0> sample_idx[n_shards, N];
	int<lower=0> G_per_shard[n_shards];
	int<lower=0> G_per_shard_idx[n_shards + 1];

	// Global properies prior model
	real lambda_mu_mu;
	real lambda_sigma;
  real<upper=0> sigma_slope;
  real sigma_intercept;
  real<lower=0>sigma_sigma;

  // Deconvolution
  //int<lower=1> C;
  int<lower=0> Q;
  int<lower=0> I;
  int<lower=1> ct_in_levels[2];
  int y[I,2]; // `read count`  S Q mapping_with_ct

  // Skip lambda imputaton
 	int<lower=0, upper=1> do_infer;
 	vector[G] lambda_log_data;
  vector[G] sigma_raw_data;

  // Efficient sum calculation
  int I1;
  int I2;
  int I1_dim[2];
  int I2_dim[2];
	int idx_1[I1];
	int idx_2[I2];

	int y_MPI_symbol_per_shard[n_shards];
	int y_MPI_G_per_shard[n_shards];
	int y_MPI_idx[n_shards,max(y_MPI_G_per_shard)];
	int y_MPI_N_per_shard[n_shards];
	int y_MPI_count[n_shards, max(y_MPI_N_per_shard)];

}
transformed data {

	int y_1_rows = I1_dim[1] * Q;
	int y_2_rows = I2_dim[1] * Q;

  vector[0] global_parameters;
  real xr[n_shards, 0];

  int<lower=0> int_MPI[n_shards, 4+(M+1)+N+N];
	int<lower=0> int_MPI_lv1[n_shards, 6 + max(y_MPI_N_per_shard)];

	/////////////////////////////////////////////
	// Shards - reference MPI
	/////////////////////////////////////////////

  for ( i in 1:n_shards ) {
	  int M_N_Gps[4];
	  M_N_Gps[1] = M;
	  M_N_Gps[2] = N;
	  M_N_Gps[3] = S;
	  M_N_Gps[4] = G_per_shard[i];

	  int_MPI[i,] = append_array(append_array(append_array(M_N_Gps, symbol_end[i]), sample_idx[i]), counts[i]);

  }

	/////////////////////////////////////////////
	// Shards - deconvolutoin MPI level 1
	/////////////////////////////////////////////


  for ( i in 1:n_shards ) {
	  int M_N_Gps[6];
	  M_N_Gps[1] = ct_in_levels[1];
	  M_N_Gps[2] = Q;
	  M_N_Gps[3] = S;
	  M_N_Gps[4] = y_MPI_symbol_per_shard[i];
	  M_N_Gps[5] = y_MPI_G_per_shard[i];
	  M_N_Gps[6] = y_MPI_N_per_shard[i];

	  int_MPI_lv1[i] = append_array(M_N_Gps, y_MPI_count[i]);

  }


}
parameters {
  // Overall properties of the data
  real<lower=0> lambda_mu; // So is compatible with logGamma prior
  //real<lower=0> lambda_sigma;
  vector[S] exposure_rate;

  // Gene-wise properties of the data
  vector[G * do_infer] lambda_log_param;
  vector[G * do_infer] sigma_raw_param;

  // Proportions
  simplex[ct_in_levels[1]] prop_1[Q]; // Root
  simplex[ct_in_levels[2]] prop_immune[Q]; // Immune cells

}
transformed parameters {
  // Sigma
  vector[G] sigma = 1.0 ./ exp(do_infer ? sigma_raw_param : sigma_raw_data) ;
	vector[G] lambda_log = do_infer ? lambda_log_param : lambda_log_data;
	vector[G] lambda = exp(lambda_log);

	// proportion of the higher level
	vector[sum(ct_in_levels) - 1] prop_2[Q] = append_vector_array(
		prop_1[, 1:3],
		multiply_by_column(prop_immune, prop_1[, 4])
	);

	vector[2*M + S] lambda_sigma_exposure_MPI[n_shards];
	vector[max(y_MPI_G_per_shard) * 2 + Q + (Q * ct_in_levels[1])] lambda_sigma_exposure_prop_MPI_lv1[n_shards];

	/////////////////////////////////////////////
	// Shards - reference MPI
	/////////////////////////////////////////////

	for( i in 1:n_shards ) {

		int size_buffer = (M*2) - (G_per_shard[i]*2) ;
	  vector[ size_buffer] buffer = rep_vector(0.0,size_buffer);

		lambda_sigma_exposure_MPI[i] =
  		append_row(
  		  append_row(
	  		    append_row(
	  		    	lambda_log[(G_per_shard_idx[i]+1):(G_per_shard_idx[i+1])],
	      		  sigma[(G_per_shard_idx[i]+1):(G_per_shard_idx[i+1])]
	      		),
      		buffer
      	),
      	exposure_rate
      );
	}

	/////////////////////////////////////////////
	// Shards - deconvolutoin MPI level 1
	/////////////////////////////////////////////

	for( i in 1:n_shards ) {

		int size_buffer = ( max(y_MPI_G_per_shard) - y_MPI_G_per_shard[i] ) * 2;
	  vector[size_buffer] buffer = rep_vector(0.0,size_buffer);

		lambda_sigma_exposure_prop_MPI_lv1[i] =
			append_row(
				append_row(
				  append_row(
					    append_row(
					    	lambda[y_MPI_idx[i, 1:y_MPI_G_per_shard[i]]],
			    		  sigma[y_MPI_idx[i, 1:y_MPI_G_per_shard[i]]]
			    		),
			  		exposure_rate[1:Q]
				  	),
				  vector_array_to_vector(prop_1)
				  ),
				buffer
			);
	}
}

model {

	//	vector[y_1_rows * 2] sum1 = sum_NB( lambda[idx_1], sigma[idx_1], I1_dim, prop_1);
	vector[y_2_rows * 2] sum2 = sum_NB( lambda[idx_2], sigma[idx_2], I2_dim, prop_2);

	target += sum( map_rect( sum_reduce , global_parameters , lambda_sigma_exposure_prop_MPI_lv1 , xr , int_MPI_lv1 ) );
	//target +=  sum_reduce (global_parameters , lambda_sigma_exposure_prop_MPI_lv1[1] , xr[1] , int_MPI_lv1[1] ) ;

	// Vecotrised sampling
	y[y_1_rows+1:size(y),1]  ~ neg_binomial_2(
		sum2[1:y_2_rows] .* exp(exposure_rate)[y[y_1_rows+1:size(y),2]] ,
	 	sum2[(y_2_rows+1):(y_2_rows*2)]
	);

  // Overall properties of the data
  lambda_mu ~ normal(lambda_mu_mu,2);

	// Exposure prior
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.001 * S);

  // Gene-wise properties of the data
  if(do_infer) lambda_log_param ~ exp_gamma_meanSd(lambda_mu,lambda_sigma);
  if(do_infer) sigma_raw_param ~ normal(sigma_slope * lambda_log_param + sigma_intercept,sigma_sigma);

	// Gene-wise properties of the data
	target += sum( map_rect( lp_reduce , global_parameters , lambda_sigma_exposure_MPI , xr , int_MPI ) );

}
// generated quantities{
//   vector[G] sigma_raw_gen;
//   vector[G] lambda_gen;
//
//   for(g in 1:G) sigma_raw_gen[g] = normal_rng(sigma_slope * lambda_log[g] + sigma_intercept,sigma_sigma);
//   for(g in 1:G) lambda_gen[g] =  exp_gamma_meanSd_rng(lambda_mu,lambda_sigma);
//
//
// }
