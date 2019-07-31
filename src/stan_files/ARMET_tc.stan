functions{

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
}
data {

	// Reference matrix inference

	int<lower=0> G;
	int<lower=0> S;
	int CL;

 	int counts_linear[CL] ;
	int G_linear[CL] ;
	int S_linear[CL] ;

	real<upper=0> sigma_slope;
  real<lower=0>sigma_sigma;

 	int<lower=0> Q;
  int<lower=1> n_nodes;
  int<lower=1> ct_in_nodes[n_nodes];
	int<lower=1> n_levels;
  int<lower=1> ct_in_levels[n_levels];

  // Deconvolution
  int<lower=0> GM;
  int GM1_linear[GM];
  int<lower=0> Y;
	int y_linear[Y];

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

}
parameters {

  // Gene-wise properties of the data
  vector[G] lambda_log;

  vector[S] exposure_rate;

	real<lower=0> lambda_mu;
  real<lower=0> lambda_sigma;
  real<upper=0> lambda_skew;

  real sigma_intercept;

  // Proportions
  simplex[ct_in_nodes[1]] prop_1[Q]; // Root
  simplex[ct_in_nodes[2]] prop_a[Q]; // Immune cells
  simplex[ct_in_nodes[3]] prop_b[Q]; // b cells
  simplex[ct_in_nodes[4]] prop_c[Q]; // granulocyte
  simplex[ct_in_nodes[5]] prop_d[Q]; // mono_derived
  simplex[ct_in_nodes[6]] prop_e[Q]; // t_cell

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

		matrix[ct_in_nodes[1], GM/ct_in_nodes[1]] mat_GM1 = exp(to_matrix(lambda_log[GM1_linear], ct_in_nodes[1], GM/ct_in_nodes[1])) ;
}

model {

	// Calculate dicunvolute
	vector[Y] lambda_log_deconvoluted =
		log(
			to_vector(
				vector_array_to_matrix(prop_1) *
				exp(to_matrix(lambda_log[GM1_linear], ct_in_nodes[1], GM/ct_in_nodes[1])) // [Q,G] dimensions
			)
		);


  // Overall properties of the data
  lambda_mu ~ normal(0,2);
	lambda_sigma ~ normal(0,2);
	lambda_skew ~ normal(0,1);
	sigma_intercept ~ student_t(8, 0, 1);

	// Exposure
	exposure_rate ~ normal(0,1);
	sum(exposure_rate) ~ normal(0, 0.001 * S);
	lambda_log ~ skew_normal(lambda_mu, lambda_sigma, lambda_skew);

	counts_linear ~ neg_binomial_2_log(lambda_log[G_linear] + exposure_rate[S_linear], 1.0 ./ exp( sigma_slope * lambda_log[G_linear] + sigma_intercept));

	// Deconvolution
	for(q in 1:Q) prop_1[q] ~ dirichlet(rep_vector(num_elements(prop_1[1]), num_elements(prop_1[1])));
	for(q in 1:Q) prop_a[q] ~ dirichlet(rep_vector(num_elements(prop_a[1]), num_elements(prop_a[1])));
	for(q in 1:Q) prop_b[q] ~ dirichlet(rep_vector(num_elements(prop_b[1]), num_elements(prop_b[1])));
	for(q in 1:Q) prop_c[q] ~ dirichlet(rep_vector(num_elements(prop_c[1]), num_elements(prop_c[1])));
	for(q in 1:Q) prop_d[q] ~ dirichlet(rep_vector(num_elements(prop_d[1]), num_elements(prop_d[1])));
	for(q in 1:Q) prop_e[q] ~ dirichlet(rep_vector(num_elements(prop_e[1]), num_elements(prop_e[1])));

	y_linear ~ neg_binomial_2_log(
		lambda_log_deconvoluted,
		1.0 ./ exp( sigma_slope * lambda_log_deconvoluted + sigma_intercept)
	);

}
