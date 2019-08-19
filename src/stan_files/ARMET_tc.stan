functions{

 vector horseshoe_get_tp(vector zb, vector[] local, real[] global, real scale_global, real c2) {
  	    int K = rows(zb);
  	    vector[K] lambda = local[1] .* sqrt(local[2]);
  	    vector[K] lambda2 = square(lambda);
  	    real tau = global[1] * sqrt(global[2]) * scale_global;
  	    vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
  	    return zb .* lambda_tilde * tau;
  	  }

  real horseshoe_get_lp(vector zb, vector[] local, real df, real[] global, real df_global, real c2, real df_slab){

    // real<lower=0> hs_df; // == 1  // If divergencies increase this
	  // real<lower=0> hs_df_global; // == 1
	  // real<lower=0> hs_df_slab; // == 4 // df of the outliers

  	vector[6] lp;

  	lp[1] = normal_lpdf(zb | 0, 1);
	  lp[2] = normal_lpdf(local[1] | 0, 1) - 101 * log(0.5);
	  lp[3] = inv_gamma_lpdf(local[2] | 0.5 * df, 0.5 * df);
	  lp[4] = normal_lpdf(global[1] | 0, 1)  - 1 * log(0.5);
	  lp[5] = inv_gamma_lpdf(global[2] | 0.5 * df_global, 0.5 * df_global);
	  lp[6] = inv_gamma_lpdf(c2 | 0.5 * df_slab, 0.5 * df_slab);

	  return(sum(lp));
  }

	matrix vector_array_to_matrix(vector[] x) {
			matrix[size(x), rows(x[1])] y;
			for (m in 1:size(x))
			  y[m] = x[m]';
			return y;
	}

	vector[] multiply_by_column(vector[] v, real[] r){
		int n_rows = num_elements(v[,1]);
		int n_cols = num_elements(v[1]);

		vector[n_cols] v_mult [n_rows];
		for(i in 1:n_cols) v_mult[,i] = to_array_1d(to_row_vector(v[,i]) .* to_row_vector(r));

		return v_mult;
	}

	vector[] append_vector_array(vector[] v1, vector[] v2){
		vector[num_elements(v1[1]) + num_elements(v2[1])] v3[num_elements(v1[,1])];

		v3[,1:num_elements(v1[1])] = v1;
		v3[,num_elements(v1[1])+1:num_elements(v3[1])] = v2;

		return v3;
	}

	int get_buffer_size(vector v, real threshold){
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

	vector[] get_real_MPI(vector v, int shards){

		int elements_per_shard[shards] = get_elements_per_shard(rows(v), shards); // Length of the returned object
		int size_MPI_obj = elements_per_shard[1]; // the first element is always the full size of the object
		vector[size_MPI_obj] v_MPI[shards] ; // Set values to -999 for the ones that are not filled

		int i = 0; // Index sweeping the vector

		for(s in 1:shards){
			v_MPI[s] = rep_vector(-999.0, size_MPI_obj);
			v_MPI[s, 1:elements_per_shard[s]] = v[ (i + 1) : i + elements_per_shard[s] ];
			i += elements_per_shard[s];
		}

		return v_MPI;
	}

	vector lp_reduce( vector global_parameters , vector local_parameters , real[] real_data , int[] int_data ) {

		real lp;
		real threshold = -999;
		int size_buffer = get_buffer_size(local_parameters, threshold);
		int size_vector = rows(local_parameters)-size_buffer;

		if(min(local_parameters[1:size_vector]) == threshold) print("ERROR! The MPI implmentation is buggy")

		// Reference / exposure rate
		lp = neg_binomial_2_log_lpmf(
			int_data[1:size_vector] |
			local_parameters[1:size_vector],
			1.0 ./ exp( global_parameters[2] * local_parameters[1:size_vector] + global_parameters[1])
		);

	 return [lp]';

	}

}
data {

	// Reference matrix inference

	int<lower=0> G;
	int<lower=0> GM;
	int<lower=0> S;
	int CL;

 	int counts_linear[CL] ;
	int G_linear[CL] ;
	int S_linear[CL] ;

	real<upper=0> sigma_slope;
	real<lower=0> sigma_sigma;
	real sigma_intercept;


 	int<lower=0> Q;
  int<lower=1> n_nodes;
  int<lower=1> ct_in_nodes[n_nodes];
	int<lower=1> n_levels;
  int<lower=1> ct_in_levels[n_levels];

  // Deconvolution
  int<lower=0> G1;
  int G1_linear[G1];
  int<lower=0> G2;
  int G2_linear[G2];
  int<lower=0> G3;
  int G3_linear[G3];

  // Observed counts
  int<lower=0> Y_1;
	int y_linear_1[Y_1];
	int y_linear_S_1[Y_1];
	int y_linear_GM_1[Y_1];
  int<lower=0> Y_2;
	int y_linear_2[Y_2];
	int y_linear_S_2[Y_2];
	int y_linear_GM_2[Y_2];
	int<lower=0> Y_3;
	int y_linear_3[Y_3];
	int y_linear_S_3[Y_3];
	int y_linear_GM_3[Y_3];

	// Lv2 tree structure parents singles
	int<lower=1> SLV2;
	int<lower=1> PLV2;
	int parents_lv2[PLV2]; // Level one parents
	int singles_lv2[SLV2]; // Level 1 leafs

	// Lv3 tree structure parents singles
	int<lower=1> SLV3;
	int<lower=1> PLV3;
	int parents_lv3[PLV3]; // Level one parents
	int singles_lv3[SLV3]; // Level 1 leafs

	// Non-centered param
	real lambda_mu_prior[2];
	real lambda_sigma_prior[2];
	real lambda_skew_prior[2];
	real sigma_intercept_prior[2];

	// MPI
	int shards;

	// Horseshoe
	real<lower=0> hs_df; // == 1  // If divergencies increase this
	real<lower=0, upper=1> par_ratio; // real<lower=0> hs_scale_global; // Ratio of the expected number of non-zero coefficients     !! KEY PARAMETER
	real<lower=0> hs_scale_slab; // == 2 // regularisation/scale of outliers                !! KEY PARAMETER
	real df_global;
	real df_slab;

}
transformed data{
	// MPI
	int data_integer[CL + Y_1 + Y_2 + Y_3] =
	append_array(
		append_array(
			append_array(	counts_linear,	y_linear_1),
			y_linear_2
		),
		y_linear_3
	);

	int y_linear[Y_1 + Y_2 + Y_3] =	append_array(append_array(	y_linear_1,	y_linear_2),	y_linear_3);

	real real_data[shards, 0] = rep_array(0.0, shards, 0);
}
parameters {

	// Global properties
	real<offset=lambda_mu_prior[1],multiplier=lambda_mu_prior[2]>lambda_mu;
  real<offset=lambda_sigma_prior[1],multiplier=lambda_sigma_prior[2]> lambda_sigma;
  real<offset=lambda_skew_prior[1],multiplier=lambda_skew_prior[2]> lambda_skew;
  //real<offset=sigma_intercept_prior[1],multiplier=sigma_intercept_prior[2]> sigma_intercept;
	real sigma_intercept_dec;
	//real<upper=0> sigma_slope_dec;

  // Local properties of the data
  vector[G] lambda_log;
  vector[G] sigma_inv_log;
  vector[S] exposure_rate;

  // Proportions
  simplex[ct_in_nodes[1]] prop_1[Q]; // Root
  simplex[ct_in_nodes[2]] prop_a[Q]; // Immune cells
  simplex[ct_in_nodes[3]] prop_b[Q]; // b cells
  simplex[ct_in_nodes[4]] prop_c[Q]; // granulocyte
  simplex[ct_in_nodes[5]] prop_d[Q]; // mono_derived
  simplex[ct_in_nodes[6]] prop_e[Q]; // t_cell

  // Error between reference and mix, to avoid divergencies
  vector<lower=0>[GM] error_ref_mix_z;

  // Horseshoe
  vector<lower=0>[GM] hs_local[2]; // local parameters for horseshoe prior
  real<lower=0> hs_global[2]; // horseshoe shrinkage parameters
  real<lower=0> hs_c2; // horseshoe shrinkage parameters

}
transformed parameters{
		// proportion of level 2
	vector[ct_in_levels[2]] prop_2[Q] =
		append_vector_array(
			prop_1[,singles_lv2],
			multiply_by_column(prop_a, prop_1[,parents_lv2[1]])
		);

			// proportion of level 3
	vector[ct_in_levels[3]] prop_3[Q] =
		append_vector_array(
			prop_2[,singles_lv3],
			append_vector_array(
				multiply_by_column(prop_b, prop_2[,parents_lv3[1]]),
				append_vector_array(
					multiply_by_column(prop_c, prop_2[,parents_lv3[2]]),
					append_vector_array(
						multiply_by_column(prop_d, prop_2[,parents_lv3[3]]),
						multiply_by_column(prop_e, prop_2[,parents_lv3[4]])
					)
				)
			)
		);

	// Horseshoe
  vector[GM] error_ref_mix = horseshoe_get_tp(error_ref_mix_z, hs_local, hs_global, par_ratio / sqrt(GM), hs_scale_slab^2 * hs_c2);


}

model {

	// Calculate convoluted
	vector[Y_1] lambda_log_deconvoluted_1 =
		log(
			to_vector(
				vector_array_to_matrix(prop_1) *
				exp(to_matrix(lambda_log[G1_linear], ct_in_levels[1], G1/ct_in_levels[1])) // [Q,G] dimensions
			)
		);

	vector[Y_2] lambda_log_deconvoluted_2 =
		log(
			to_vector(
				vector_array_to_matrix(prop_2) *
				exp(to_matrix(lambda_log[G2_linear], ct_in_levels[2], G2/ct_in_levels[2])) // [Q,G] dimensions
			)
		);

	vector[Y_3] lambda_log_deconvoluted_3 =
		log(
			to_vector(
				vector_array_to_matrix(prop_3) *
				exp(to_matrix(lambda_log[G3_linear], ct_in_levels[3], G3/ct_in_levels[3])) // [Q,G] dimensions
			)
		);

		vector[Y_1 + Y_2 + Y_3] lambda_log_deconvoluted =
			append_row(
				append_row(
					lambda_log_deconvoluted_1 + exposure_rate[y_linear_S_1],
					lambda_log_deconvoluted_2 + exposure_rate[y_linear_S_2]
				),
				lambda_log_deconvoluted_3 + exposure_rate[y_linear_S_3]
			);

	vector[Y_1 + Y_2 + Y_3] sigma_inv_log_deconvoluted = error_ref_mix[append_array(y_linear_GM_1, append_array(y_linear_GM_2, y_linear_GM_3))];

	// vector[CL + Y_1 + Y_2 + Y_3] local_parameters =
	// 	append_row(
	// 		append_row(
	// 			append_row(
	// 				lambda_log[G_linear] + exposure_rate[S_linear],
	// 				lambda_log_deconvoluted_1 + exposure_rate[y_linear_S_1]
	// 			),
	// 			lambda_log_deconvoluted_2 + exposure_rate[y_linear_S_2]
	// 		),
	// 		lambda_log_deconvoluted_3 + exposure_rate[y_linear_S_3]
	// 	);

  // Overall properties of the data
  lambda_mu ~ normal(lambda_mu_prior[1],lambda_mu_prior[2]);
	lambda_sigma ~ normal(lambda_sigma_prior[1],lambda_sigma_prior[2]);
	lambda_skew ~ normal(lambda_skew_prior[1],lambda_skew_prior[2]);
	//sigma_intercept ~ normal(sigma_intercept_prior[1], sigma_intercept_prior[2]);
	//sigma_sigma ~ normal(0,1);

	// Exposure
	exposure_rate ~ normal(0,1);
	sum(exposure_rate) ~ normal(0, 0.001 * S);

	// Means overdispersion reference
	lambda_log ~ skew_normal(lambda_mu, exp(lambda_sigma), lambda_skew);
	sigma_inv_log ~ normal(sigma_slope * lambda_log + sigma_intercept, sigma_sigma);

	// Deconvolution
	for(q in 1:Q) prop_1[q] ~ dirichlet(rep_vector(num_elements(prop_1[1]), num_elements(prop_1[1])));
	for(q in 1:Q) prop_a[q] ~ dirichlet(rep_vector(num_elements(prop_a[1]), num_elements(prop_a[1])));
	for(q in 1:Q) prop_b[q] ~ dirichlet(rep_vector(num_elements(prop_b[1]), num_elements(prop_b[1])));
	for(q in 1:Q) prop_c[q] ~ dirichlet(rep_vector(num_elements(prop_c[1]), num_elements(prop_c[1])));
	for(q in 1:Q) prop_d[q] ~ dirichlet(rep_vector(num_elements(prop_d[1]), num_elements(prop_d[1])));
	for(q in 1:Q) prop_e[q] ~ dirichlet(rep_vector(num_elements(prop_e[1]), num_elements(prop_e[1])));

	// Unexplanable difference between reference and mix
	// error_ref_mix ~ exponential(1);
	// Horseshoe
  target += horseshoe_get_lp(error_ref_mix_z, hs_local, hs_df, hs_global, df_global, hs_c2, df_slab);


	// target += sum(map_rect(
	// 	lp_reduce ,
	// 	[sigma_intercept, sigma_slope]', // global parameters
	// 	get_real_MPI(local_parameters, shards),
	// 	real_data, // real data
	// 	get_int_MPI(data_integer, shards)
	// ));

	// Deconvolution
	sigma_intercept_dec ~ student_t(3, 0, 2);
	//sigma_slope_dec ~ normal(0,1);

	y_linear ~ neg_binomial_2_log(
		lambda_log_deconvoluted,
		1.0 ./ exp( sigma_inv_log_deconvoluted + ( sigma_slope * lambda_log_deconvoluted + sigma_intercept_dec ) )
	);

	// Reference
	counts_linear ~ neg_binomial_2_log(
		lambda_log[G_linear] + exposure_rate[S_linear],
		1.0 ./ exp( sigma_inv_log[G_linear] )
	);

}
