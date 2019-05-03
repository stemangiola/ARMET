functions{

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

}
transformed data {

  vector[0] global_parameters;
  real xr[n_shards, 0];

  int<lower=0> int_MPI[n_shards, 4+(M+1)+N+N];

  // Shards - MPI
  for ( i in 1:n_shards ) {
  int M_N_Gps[4];
  M_N_Gps[1] = M;
  M_N_Gps[2] = N;
  M_N_Gps[3] = S;
  M_N_Gps[4] = G_per_shard[i];

  int_MPI[i,] = append_array(append_array(append_array(M_N_Gps, symbol_end[i]), sample_idx[i]), counts[i]);

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

	// proportion of the higher level
	vector[sum(ct_in_levels) - 1] prop_2[Q] =
		append_vector_array(	prop_1[, 1:3], multiply_by_column(prop_immune, prop_1[, 4]) );

	// Shards - MPI
	vector[2*M + S] lambda_sigma_exposure_MPI[n_shards];
	for( i in 1:(n_shards) ) {

	  vector[ (M*2) - (G_per_shard[i]*2) ] buffer = rep_vector(0.0,(M*2) - (G_per_shard[i]*2));

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
}

model {

	// Deconvolution
	// prop_1 ~ dirichlet(rep_vector(1*ct_in_levels[1], ct_in_levels[1]));
	// prop_immune ~ dirichlet(rep_vector(1*ct_in_levels[2], ct_in_levels[2]));

	vector[G] lambda = exp(lambda_log);



	matrix[I1_dim[1], I1_dim[2]] lambda_mat_1 = to_matrix(lambda[idx_1], I1_dim[1], I1_dim[2]); // Bug prone
	matrix[rows(prop_1[1]),Q] prop_mat_1 = vector_array_to_transposed_matrix(prop_1);
	matrix[I1_dim[1] , Q] lambda_sum_1 = lambda_mat_1 * prop_mat_1;
	//matrix[I1_dim[1], I1_dim[2]] sigma_mat_1 = to_matrix(sigma[idx_1], I1_dim[1], I1_dim[2]);
	matrix[I1_dim[1], Q] sigma_sum_1 =
		square(lambda_sum_1) ./
		((
			square(lambda_mat_1) ./
			to_matrix(sigma[idx_1], I1_dim[1], I1_dim[2])
		) * square(prop_mat_1)) ;


	matrix[I2_dim[1], I2_dim[2]] lambda_mat_2 = to_matrix(lambda[idx_2], I2_dim[1], I2_dim[2]); // Bug prone
	matrix[rows(prop_2[1]), Q ] prop_mat_2 = vector_array_to_transposed_matrix(prop_2);
	matrix[I2_dim[1] , Q] lambda_sum_2 = lambda_mat_2 * prop_mat_2;
	//matrix[I2_dim[1], I2_dim[2]] sigma_mat_2 = to_matrix(sigma[idx_2], I2_dim[1], I2_dim[2]);
	matrix[I2_dim[1], Q] sigma_sum_2 =
		square(lambda_sum_2) ./
		((
			square(lambda_mat_2) ./
			to_matrix(sigma[idx_2], I2_dim[1], I2_dim[2])
		) * square(prop_mat_2));


	// Vecotrised sampling
	y[,1]  ~ neg_binomial_2(
		append_row( to_vector(lambda_sum_1), to_vector(lambda_sum_2)) .* exp(exposure_rate)[y[,2]] ,
		append_row( to_vector(sigma_sum_1), to_vector(sigma_sum_2))
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
