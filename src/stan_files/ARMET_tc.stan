functions{

  	real exp_gamma_meanSd_lpdf(vector x_log, real m_log, real s){

      // This function is the  probability of the log gamma function
      // in case you have data that is aleady in log form

      // real v = square(s);
      // real a = square(m) / v;
      // real b = m / v;

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

	vector get_sum_NB(vector mu, vector phi){

			// Conversion to gama
			vector[rows(mu)] shape = (mu .* phi)./(mu+phi);
			vector[rows(mu)] rate = phi ./ (mu+phi);

			// Calculating gamma
			real sum_shape_on_rate = sum( shape ./ rate );
			real shape_sum = sum_shape_on_rate^2/sum(shape ./ (rate .* rate));
			real rate_sum = 1/(sum_shape_on_rate / shape_sum);

			// Conversion to NB

			// mu_sum
			real mu_sum = shape_sum / rate_sum; // or sum(mu);
			// sigma_sum
			real sigma_sum = mu_sum^2 / ( (shape_sum/rate_sum^2) - mu_sum );

			vector[2] mu_sigma;
			mu_sigma[1] = mu_sum;
			mu_sigma[2] = sigma_sum;

			return mu_sigma;
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
  int<lower=1> ct_in_levels[1];
  int y[I, 3 + max(ct_in_levels)]; // `read count`  S Q mapping_with_ct

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
  vector[G] lambda_log;
  vector[G] sigma_raw;

  // Proportions
  simplex[ct_in_levels[1]] prop[Q];

}
transformed parameters {
  // Sigma
  vector[G] sigma = 1.0 ./ exp(sigma_raw);

// Shards - MPI
	vector[2*M] lambda_sigma_MPI[n_shards];
	for( i in 1:(n_shards) ) {

	  vector[ (M*2) - (G_per_shard[i]*2) ] buffer = rep_vector(0.0,(M*2) - (G_per_shard[i]*2));

		lambda_sigma_MPI[i] =
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

	vector[2] y_sum[I];

  // Overall properties of the data
  lambda_mu ~ normal(lambda_mu_mu,2);
  //lambda_sigma ~ normal(0,2);

  //sigma_raw ~ normal(0,1);
  exposure_rate ~ normal(0,1);
  sum(exposure_rate) ~ normal(0, 0.001 * S);

  // Gene-wise properties of the data
  // lambda_log ~ normal_or_gammaLog(lambda_mu, lambda_sigma, is_prior_asymetric);
  lambda_log ~ exp_gamma_meanSd(lambda_mu,lambda_sigma);

  sigma_raw ~ normal(sigma_slope * lambda_log + sigma_intercept,sigma_sigma);

	// Gene-wise properties of the data
	target += sum( map_rect( lp_reduce , global_parameters , lambda_sigma_MPI , xr , int_MPI ) );

	// Deconvolution

	for(i in 1:I) y_sum[i] =
		get_sum_NB(
			exp(lambda_log[y[i,4 : 3 + ct_in_levels[1]]]) .* prop[y[i,3]],
			sigma[y[i,4 : 3 + ct_in_levels[1]]]
		);

	y[,1]  ~ neg_binomial_2_log(to_vector( log( y_sum[,1]) ) + exposure_rate[y[,2]] , y_sum[,2]);

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