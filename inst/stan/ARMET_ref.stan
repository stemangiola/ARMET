functions{

int[] rep_int(int x, int n){
	int out_int[n];

	for(i in 1:n)	out_int[i] = x;

	return out_int;
}

// takes as input an array of the form {a,b,c,d} where a,b,c,d are
// real 1d arrays. These are then returned as concatenated 1d array
real[] concatenate_real_array(real[,] elems);
real[] concatenate_real_array(real[,] elems) {
  int num_elems = size(elems);

  if (num_elems == 1)
    return(elems[1]);

  if (num_elems == 2)
    return(to_array_1d(append_row(to_vector(elems[1]), to_vector(elems[2]))));

  // else
  return(to_array_1d(append_row(to_vector(elems[1]), to_vector( concatenate_real_array(elems[2:num_elems]) ))));
}

// takes as input an array of the form {a,b,c,d} where a,b,c,d are
// 1d vectors. These are then returned as concatenated 1d vector
vector concatenate_vector_array(vector[] elems);
vector concatenate_vector_array(vector[] elems) {
  int num_elems = size(elems);

  if (num_elems == 1) return(elems[1]);

  if (num_elems == 2) return(append_row(elems[1], elems[2]));

  // else
 	return(append_row(elems[1],  concatenate_vector_array(elems[2:num_elems]) ));
}

int[] append_int(int[] i1, int[] i2){
	int i3[size(i1)+size(i2)];

	i3[1:size(i1)] = i1;
	i3[size(i1)+1:size(i1)+size(i2)] = i2;

	return i3;
}

// takes as input an array of the form {a,b,c,d} where a,b,c,d are
// 1d int[. These are then returned as concatenated 1d integer
int[] concatenate_int_array(int[,] elems);
int[] concatenate_int_array(int[,] elems) {
  int num_elems = size(elems);

  if (num_elems == 1) return(elems[1]);

  if (num_elems == 2) return(append_int(elems[1], elems[2]));

  // else
 	return(append_int(elems[1],  concatenate_int_array(elems[2:num_elems]) ));
}

int get_real_buffer_size(vector v, real threshold){
	// This function finds how may fake indexes -1 there are in a vector, added for map_rect needs

	real i = threshold; // Value of the index
	int n = 0; // Length of the buffer
	int s = rows(v); // Size of the whole vector

	while(i == threshold){
		i = v[s-n];
		if(i==threshold) n += 1;
	}

	return n;
}

int[] get_elements_per_shard(int lenth_v, int shards){

	// Returned integer(max_size, last_element_size)
	int tentative_size = lenth_v / shards;
	int tentative_remaining = lenth_v - (tentative_size * shards);
	int elements_per_shard = tentative_remaining > 0 ? tentative_size + 1 : tentative_size;
	int remaining =  (elements_per_shard * shards) - lenth_v;

	int length_obj[shards];

	for(s in 1:shards) {
		length_obj[s] =
			s != shards ?
			elements_per_shard :
			elements_per_shard - remaining;  // Actual elements in it for last object
	}

 	return length_obj;

}

int[,] get_int_MPI(int[] v, int shards){
  // Simple MPI for int vector

	int elements_per_shard[shards] = get_elements_per_shard(size(v), shards); // Length of the returned object
	int size_MPI_obj = elements_per_shard[1]; // the first element is always the full size of the object
	int v_MPI[shards,size_MPI_obj] = rep_array(-1, shards,size_MPI_obj); // Set values to -1 for the ones that are not filled

	int i = 0; // Index sweeping the vector

	for(s in 1:shards){
		v_MPI[s, 1:elements_per_shard[s]] = v[ (i + 1) : i + elements_per_shard[s] ];
		i += elements_per_shard[s];
	}

	return v_MPI;
}

	vector[] get_mu_sigma_vector_MPI(vector mus, vector sigmas, int shards){

		int elements_per_shard[shards] = get_elements_per_shard(rows(mus), shards); // Length of the returned object
		int size_MPI_obj = elements_per_shard[1]; // the first element is always the full size of the object
		vector[size_MPI_obj * 2] v_MPI[shards] ; // Set values to -999 for the ones that are not filled

		int i = 0; // Index sweeping the vector

		for(s in 1:shards){

			// If last shard fill in
			if(s == shards) v_MPI[s] = rep_vector(-999.0, size_MPI_obj * 2);

			v_MPI[s, 1:elements_per_shard[s]] = mus[ (i + 1) : i + elements_per_shard[s] ];
			v_MPI[s, (elements_per_shard[s]+1):(elements_per_shard[s]+elements_per_shard[s])] = sigmas[ (i + 1) : i + elements_per_shard[s] ];

			i += elements_per_shard[s];
		}

		return v_MPI;
	}

	vector lp_reduce_simple( vector global_parameters , vector mus_sigmas , real[] real_data , int[] int_data ) {

		real lp;
		real threshold = -999;
		int size_buffer = get_real_buffer_size(mus_sigmas, threshold);
		int size_vector = (rows(mus_sigmas)-size_buffer)/2;

		if(min(mus_sigmas[1:(size_vector*2)]) == threshold) print("ERROR! The MPI implmentation is buggy")

		// Reference / exposure rate
		lp = neg_binomial_2_log_lpmf(
			int_data[1:size_vector] |
			mus_sigmas[1:size_vector],
			1.0 ./ exp( mus_sigmas[size_vector+1:size_vector+size_vector] )
		);

	 return [lp]';

	}

	vector[] which(int x, vector[] a, vector[] b, vector[] c, vector[] d){
		if(x == 1) return(a);
		if(x == 2) return(b);
		if(x == 3) return(c);
		else return(d);
	}


}
data {
	// shards
	int<lower=1> shards;

	// Reference matrix inference
	int<lower=0> G;
	int<lower=0> GM;
	int<lower=0> S;
	int CL; // counts linear size

	// reference counts
 	int<lower=0> counts_linear[CL] ;
	int G_to_counts_linear[CL] ;
	int S_linear[CL] ;

	// MPI
	int<lower=0> size_counts_idx_lv_MPI[max(shards, 1)];
	int counts_idx_lv_MPI[shards, max(size_counts_idx_lv_MPI)];
	int size_counts_G_lv_MPI[max(shards, 1)];
	int size_counts_S_lv_MPI[max(shards, 1)];
	int counts_G_lv_MPI[shards,max(size_counts_G_lv_MPI)];
	int counts_S_lv_MPI[shards,max(size_counts_S_lv_MPI)];
	int size_G_linear_MPI[max(shards, 1)];
	int G_linear_MPI[shards,max(size_G_linear_MPI)];
	int size_counts_G_lv_MPI_non_redundant[max(shards, 1)];
	int counts_G_lv_MPI_non_redundant[shards, max(size_counts_G_lv_MPI_non_redundant)];
	int counts_G_lv_MPI_non_redundant_reps[shards, max(size_counts_G_lv_MPI_non_redundant)];

	// Non-centered param
	real lambda_mu_prior[2];
	real lambda_sigma_prior[2];
	real lambda_skew_prior[2];
	real sigma_intercept_prior[2];

}
transformed data{

	real real_data[shards, 0];
	real real_data2[shards, 0];

}
parameters {

	// Global properties
	real<offset=lambda_mu_prior[1],multiplier=lambda_mu_prior[2]>lambda_mu;
  real<offset=lambda_sigma_prior[1],multiplier=lambda_sigma_prior[2]> lambda_sigma;
  real<offset=lambda_skew_prior[1],multiplier=lambda_skew_prior[2]> lambda_skew;

	// Sigma
	real<upper=0> sigma_slope;
	real<lower=0> sigma_sigma;
  real<offset=sigma_intercept_prior[1],multiplier=sigma_intercept_prior[2]> sigma_intercept;

  // Local properties of the data
  vector[G] lambda_log;
  vector[G] sigma_inv_log;
  vector[S-1] exposure_rate_minus_1;


}
transformed parameters{

	vector[S] exposure_rate = append_row(exposure_rate_minus_1, -sum(exposure_rate_minus_1));

}
model {

  // Overall properties of the data
  lambda_mu ~ normal(lambda_mu_prior[1],lambda_mu_prior[2]);
	lambda_sigma ~ normal(lambda_sigma_prior[1],lambda_sigma_prior[2]);
	lambda_skew ~ normal(lambda_skew_prior[1],lambda_skew_prior[2]);

  sigma_intercept ~ normal(0,2);
  sigma_slope ~ normal(0,2);
  sigma_sigma ~ normal(0,2);

	// Exposure
	exposure_rate_minus_1 ~ normal(0,1);

	// Means overdispersion reference
	lambda_log ~ skew_normal(lambda_mu, exp(lambda_sigma), lambda_skew);
	sigma_inv_log ~ normal(sigma_slope * lambda_log + sigma_intercept, sigma_sigma);

	target += sum(map_rect(
		lp_reduce_simple ,
		[sigma_intercept, sigma_slope]', // global parameters
		get_mu_sigma_vector_MPI(
			lambda_log[G_to_counts_linear] + exposure_rate[S_linear],
			sigma_inv_log[G_to_counts_linear],
			shards
		),
		real_data,
		get_int_MPI( counts_linear, shards)
	));

}
