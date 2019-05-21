functions{

		vector vector_array_to_vector(vector[] x) {
					// This operation is column major
				  vector[rows(x[1]) * size(x)] y;
				  for (m in 1:rows(x[1]))
				    y[(m-1)*size(x)+1:(m*size(x))] = to_vector(x[,m]);
				  return y;
		}

		matrix vector_array_to_transposed_matrix(vector[] x) {
		  matrix[rows(x[1]), size(x)] y;
		  for (m in 1:size(x))
		    y[,m] = x[m];
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

	vector sum_NB_MPI(matrix lambda_mat, matrix sigma_mat, matrix prop_mat){

			matrix[rows(prop_mat), cols(lambda_mat)] lambda_sum = prop_mat * lambda_mat; //Q rows, G columns
			matrix[rows(prop_mat), cols(sigma_mat)] sigma_sum =                          //Q rows, G columns
				square(lambda_sum) ./
				(
					square(prop_mat) *
					(
						square(lambda_mat) ./
						sigma_mat
					)
				) ;

				// The vectorisation is G1-Q1, G1-Q2, G1-Q3 etc..
				return(append_row( to_vector(lambda_sum), to_vector(sigma_sum)));
	}

	vector[] get_reference_parameters_MPI(int n_shards, int M, int[] G_per_shard, int[,] G_ind, vector lambda_log, vector sigma, vector exposure_rate){

		vector[2*M + cols(exposure_rate)] lambda_sigma_exposure_MPI[n_shards];

		for( i in 1:n_shards ) {

			int size_buffer = (M*2) - (G_per_shard[i]*2) ;
		  vector[ size_buffer] buffer = rep_vector(0.0,size_buffer);

			lambda_sigma_exposure_MPI[i] =
	  		append_row(
	  		  append_row(
		  		    append_row(
		  		    	lambda_log[G_ind[i, 1:G_per_shard[i]]],
		      		  sigma[G_ind[i, 1:G_per_shard[i]]]
		      		),
	      		buffer
	      	),
	      	exposure_rate
	      );
		}

		return(lambda_sigma_exposure_MPI);
	}

	vector[] get_deconvolution_parameters_MPI(
		int n_shards, int[] y_MPI_G_per_shard, int[,] y_MPI_idx, vector lambda, vector sigma, vector exposure_rate, vector sigma_raw_global, int Q, vector[] prop){

		vector[max(y_MPI_G_per_shard) * 2 + Q + Q + (Q * cols(prop[1]))] lambda_sigma_exposure_prop_MPI[n_shards];

		for( i in 1:n_shards ) {

			int size_buffer = ( max(y_MPI_G_per_shard) - y_MPI_G_per_shard[i] ) * 2;
		  vector[size_buffer] buffer = rep_vector(0.0,size_buffer);

			lambda_sigma_exposure_prop_MPI[i] =
				append_row(
					append_row(
					  append_row(
					  	append_row(
						    append_row(
						    	lambda[y_MPI_idx[i, 1:y_MPI_G_per_shard[i]]],
				    		  sigma[y_MPI_idx[i, 1:y_MPI_G_per_shard[i]]]
				    		),
					  		exposure_rate[1:Q]
						  ),
						  sigma_raw_global
						),
					  vector_array_to_vector(prop)
					),
					buffer
				);
	}

	return(lambda_sigma_exposure_prop_MPI);

}
	real reference_reduce( vector global_parameters , vector local_parameters , real[] xr , int[] xi ) {

		// Data
	 	int M = xi[1];
	 	int N = xi[2];
	 	int S = xi[3];
	 	int G_per_shard = xi[4];
	 	int symbol_end[M+1] = xi[(4+1):(4+1+M)];
	 	int sample_idx[N] = xi[(4+1+M+1):(4+1+M+1+N-1)];
	 	int counts[N] = xi[(4+1+M+1+N):size(xi)];

		// Parameters
	 	vector[G_per_shard] lambda_MPI = local_parameters[1:G_per_shard];
	 	vector[G_per_shard] sigma_MPI = local_parameters[(G_per_shard+1):(G_per_shard*2)];
	 	vector[S] exposure_rate = local_parameters[((M*2)+1):rows(local_parameters)];

		// Vectorise lpmf
		vector[symbol_end[G_per_shard+1]] lambda_MPI_c;
		vector[symbol_end[G_per_shard+1]] sigma_MPI_c;
		for(g in 1:G_per_shard){
			int how_many = symbol_end[g+1] - (symbol_end[g]);
			lambda_MPI_c[(symbol_end[g]+1):symbol_end[g+1]] = rep_vector(lambda_MPI[g], how_many);
			sigma_MPI_c [(symbol_end[g]+1):symbol_end[g+1]] = rep_vector(sigma_MPI[g],  how_many);
		}

		// Return
    return (neg_binomial_2_log_lpmf(
    	counts[1:symbol_end[G_per_shard+1]] |
    	exposure_rate[sample_idx[1:symbol_end[G_per_shard+1]]] +
    	lambda_MPI_c,
    	sigma_MPI_c
    ) * 10);

  }

  real sum_reduce( vector global_parameters , vector local_parameters , real[] xr , int[] xi ) {

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
	 	vector[Q] sigma_raw_global = local_parameters[(((y_MPI_G_per_shard*2)) + Q+1):(((y_MPI_G_per_shard*2)) + Q+Q)];

	 	vector[Q * ct_in_levels] prop = local_parameters[((y_MPI_G_per_shard*2)) + Q+Q+1	:	(((y_MPI_G_per_shard*2)) + Q+Q) + (Q * ct_in_levels)];


		vector[y_MPI_N_per_shard * 2] my_sum = sum_NB_MPI(
			to_matrix( lambda_MPI, ct_in_levels, y_MPI_symbol_per_shard), // ct rows, G columns
			to_matrix( sigma_MPI,  ct_in_levels, y_MPI_symbol_per_shard), // ct rows, G columns
			to_matrix( prop, Q, ct_in_levels)
		);

		// Vecotrised sampling, all vectors should be G1-Q1, G1-Q2, G1-Q3
		return (neg_binomial_2_lpmf(
				counts |
				my_sum[1:y_MPI_N_per_shard] .*  to_vector(rep_matrix((exp(exposure_rate)), y_MPI_symbol_per_shard)),
				my_sum[(y_MPI_N_per_shard+1):(y_MPI_N_per_shard*2)] .*  to_vector(rep_matrix((1 ./ exp(sigma_raw_global)), y_MPI_symbol_per_shard))
			));
}

	vector lp_reduce( vector global_parameters , vector local_parameters , real[] xr , int[] xi ) {

		int dim_data[3] = xi[1:3];
		int dim_param[3] = xi[4:6];

		vector[3] lp;

		lp[1] = reference_reduce(global_parameters , local_parameters[1:dim_param[1]] , xr , xi[7:(6+dim_data[1])] );
		lp[2] = sum_reduce(
			global_parameters ,
			local_parameters[(dim_param[1]+1):(dim_param[1] + dim_param[2])] ,
			xr ,
			xi[(6+dim_data[1]+1):(6+dim_data[1] + dim_data[2])]
		);
		lp[3] = sum_reduce(
			global_parameters ,
			local_parameters[(dim_param[1] + dim_param[2] +1): (dim_param[1] + dim_param[2] + dim_param[3])] ,
			xr ,
			xi[(6+dim_data[1] + dim_data[2] + 1):(6+dim_data[1] + dim_data[2] + dim_data[3])]
		);

	 return [sum(lp)]';

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
	int<lower=0> G_ind[n_shards, M];
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

  int<lower=0> counts_package[n_shards, 4+(M+1)+N+N];

	// Level 1 MPI
	int y_MPI_symbol_per_shard_lv1[n_shards];
	int y_MPI_G_per_shard_lv1[n_shards];
	int y_MPI_idx_lv1[n_shards,max(y_MPI_G_per_shard_lv1)];
	int y_MPI_N_per_shard_lv1[n_shards];
	int y_MPI_count_lv1[n_shards, max(y_MPI_N_per_shard_lv1)];
	int<lower=0> lev1_package[n_shards, 6 + max(y_MPI_N_per_shard_lv1)];

	// Level 2 MPI
	int y_MPI_symbol_per_shard_lv2[n_shards];
	int y_MPI_G_per_shard_lv2[n_shards];
	int y_MPI_idx_lv2[n_shards,max(y_MPI_G_per_shard_lv2)];
	int y_MPI_N_per_shard_lv2[n_shards];
	int y_MPI_count_lv2[n_shards, max(y_MPI_N_per_shard_lv2)];
	int<lower=0> lev2_package[n_shards, 6 + max(y_MPI_N_per_shard_lv2)];

	int data_package[n_shards, 3 + 3 + (4+(M+1)+N+N) + (6 + max(y_MPI_N_per_shard_lv1)) + (6 + max(y_MPI_N_per_shard_lv2))];
}
transformed data {

	int y_1_rows = I1_dim[1] * Q;
	int y_2_rows = I2_dim[1] * Q;

  vector[0] global_parameters;
  real xr[n_shards, 0];

}
parameters {
  // Overall properties of the data
  real<lower=0> lambda_mu; // So is compatible with logGamma prior
  //real<lower=0> lambda_sigma;
  vector[S] exposure_rate;

  // Gene-wise properties of the data
  vector[G * do_infer] lambda_log_param;
  vector[G * do_infer] sigma_raw_param;
  vector[Q] sigma_raw_global;

  // Proportions
  simplex[ct_in_levels[1]] prop_1[Q]; // Root
  simplex[ct_in_levels[2]] prop_immune[Q]; // Immune cells


}
transformed parameters {
  // Sigma
  vector[G] sigma = 1.0 ./ exp(do_infer ? sigma_raw_param : sigma_raw_data ) ;
	vector[G] lambda_log = do_infer ? lambda_log_param : lambda_log_data;

	// proportion of the higher level
	vector[sum(ct_in_levels) - 1] prop_2[Q] =
		append_vector_array(	prop_1[, 1:3], multiply_by_column(prop_immune, prop_1[, 4]) );

	// Deconvolution
	vector[G] lambda = exp(lambda_log);

}

model {

	target += sum(map_rect(
		lp_reduce ,
		global_parameters ,
		append_vector_array(
			append_vector_array(
				get_reference_parameters_MPI(	n_shards,	M, G_per_shard,	G_ind,	lambda_log,	sigma,	exposure_rate),
				get_deconvolution_parameters_MPI(n_shards, y_MPI_G_per_shard_lv1, y_MPI_idx_lv1, lambda, sigma, exposure_rate, sigma_raw_global, Q, prop_1)
			),
		get_deconvolution_parameters_MPI(n_shards, y_MPI_G_per_shard_lv2, y_MPI_idx_lv2, lambda, sigma, exposure_rate, sigma_raw_global, Q, prop_2)
		),
		xr,
		data_package
	));




	// target += sum( map_rect( sum_reduce , global_parameters , lambda_sigma_exposure_prop_MPI_lv1 , xr , lev1_package ) );
	// target += sum( map_rect( sum_reduce , global_parameters , lambda_sigma_exposure_prop_MPI_lv2 , xr , lev2_package ) );
	//

  // Overall properties of the data
  lambda_mu ~ normal(lambda_mu_mu,2);
	sigma_raw_global ~ normal(0,1);

	// Exposure prior
  exposure_rate ~ normal(0,1);
  if(do_infer) sum(exposure_rate) ~ normal(0, 0.001 * S);

  // Gene-wise properties of the data
  if(do_infer) lambda_log_param ~ exp_gamma_meanSd(lambda_mu,lambda_sigma);
  if(do_infer) sigma_raw_param ~ normal(sigma_slope * lambda_log_param + sigma_intercept,sigma_sigma);

	// // Gene-wise properties of the data
	// target += sum( map_rect( lp_reduce , global_parameters , lambda_sigma_exposure_MPI , xr , counts_package ) );



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
