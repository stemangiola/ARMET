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

int get_int_buffer_size(int[] v, int threshold){
	// This function finds how may fake indexes -1 there are in a vector, added for map_rect needs

	int i = threshold; // Value of the index
	int n = 0; // Length of the buffer
	int s = num_elements(v); // Size of the whole vector

	while(i == threshold){
		i = v[s-n];
		if(i==threshold) n += 1;
	}

	return n;
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

vector clear_real_from_buffer(vector v, real buffer_value){
	return v[1:(rows(v)- get_real_buffer_size(v, buffer_value))];
}

int[] clear_int_from_buffer(int[] v, int buffer_value){
	return v[1:(num_elements(v) - get_int_buffer_size(v, buffer_value))];
}

vector rep_vector_by_array(vector v, int[] reps){
	vector[sum(reps)] v_rep;
	int i = 0;

	for(n in 1:num_elements(reps)){
		v_rep[i+1:i+reps[n]] = rep_vector(v[n], reps[n]);
		i += reps[n];
	}

	return(v_rep);
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

vector[] get_vector_MPI(vector v, int shards){

  // Simple MPI for real vector
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

real[,] get_real_MPI(vector v, int shards){

	int elements_per_shard[shards] = get_elements_per_shard(rows(v), shards); // Length of the returned object
	int size_MPI_obj = elements_per_shard[1]; // the first element is always the full size of the object
	real v_MPI[shards,size_MPI_obj] ; // Set values to -999 for the ones that are not filled

	int i = 0; // Index sweeping the vector

	for(s in 1:shards){
		v_MPI[s,] = rep_array(-999.0, size_MPI_obj);
		v_MPI[s, 1:elements_per_shard[s]] = to_array_1d(v[ (i + 1) : i + elements_per_shard[s] ]);
		i += elements_per_shard[s];
	}

	return v_MPI;
}

int[,] append_int_MPI_arrays(int[,] lv1, int[,] lv2, int[,] lv3, int[,] lv4){

	# This is for the BUG that dim(int empty_variable[0,0]) is not {0,0} but {0}
	int dim_1[2] = size(dims(lv1)) == 2 ? dims(lv1) : rep_array(0,2);
	int dim_2[2] = size(dims(lv2)) == 2 ? dims(lv2) : rep_array(0,2);
	int dim_3[2] = size(dims(lv3)) == 2 ? dims(lv3) : rep_array(0,2);
	int dim_4[2] = size(dims(lv4)) == 2 ? dims(lv4) : rep_array(0,2);

	int max_cols = max({dim_1[2], dim_2[2],dim_3[2],dim_4[2]});
	int tot_rows = dim_1[1] + dim_2[1] + dim_3[1] + dim_4[1];


	int merged[tot_rows, max_cols];

	int i = 0;


	for(r in 1:dim_1[1]) merged[i+r] = append_array(lv1[r], rep_int(-999, max_cols-dim_1[2]));
	i += dim_1[1];

	if(dim_2[1] > 0) {
		for(r in 1:dim_2[1]) merged[i+r] = append_array(lv2[r], rep_int(-999, max_cols-dim_2[2]));
		i += dim_2[1];
	}

	if(dim_3[1] > 0) {
		for(r in 1:dim_3[1]) merged[i+r] = append_array(lv3[r], rep_int(-999, max_cols-dim_3[2]));
		i += dim_3[1];
	}

	if(dim_4[1] > 0) {
		for(r in 1:dim_4[1]) merged[i+r] = append_array(lv4[r], rep_int(-999, max_cols-dim_4[2]));
	}

	return(merged);

}

vector[] append_vector_MPI_arrays(vector[] lv1, vector[] lv2, vector[] lv3, vector[] lv4){

  // Function for appending many vectors[] inserting the buffer
	int dim_1[2];
	int dim_2[2];
	int dim_3[2];
	int dim_4[2];

	int max_cols = max({
		num_elements(lv1) > 0 ? num_elements(lv1[1]) : 0,
		num_elements(lv2) > 0 ? num_elements(lv2[1]) : 0,
		num_elements(lv3) > 0 ? num_elements(lv3[1]) : 0,
		num_elements(lv4) > 0 ? num_elements(lv4[1]) : 0
	});
	int tot_rows =
	( num_elements(lv1) > 0 ? num_elements(lv1[,1]) : 0) +
	(	num_elements(lv2) > 0 ? num_elements(lv2[,1]) : 0) +
	(	num_elements(lv3) > 0 ? num_elements(lv3[,1]) : 0) +
	(	num_elements(lv4) > 0 ? num_elements(lv4[,1]) : 0);

	vector[max_cols] merged[tot_rows];
	int i = 0;

	if(num_elements(lv1) > 0) dim_1 =  {num_elements(lv1[,1]), num_elements(lv1[1])};
	else dim_1 = rep_array(0,2);

	if(num_elements(lv2) > 0) dim_2 =  {num_elements(lv2[,1]), num_elements(lv2[1])};
	else dim_2 = rep_array(0,2);

	if(num_elements(lv3) > 0) dim_3 =  {num_elements(lv3[,1]), num_elements(lv3[1])};
	else dim_3 = rep_array(0,2);

	if(num_elements(lv4) > 0) dim_4 =  {num_elements(lv4[,1]), num_elements(lv4[1])};
	else dim_4 = rep_array(0,2);

	for(r in 1:dim_1[1]) merged[i+r] = append_row(lv1[r], rep_vector(-999, max_cols-dim_1[2]));
	i += dim_1[1];

if(dim_2[1] > 0) {
	for(r in 1:dim_2[1]) merged[i+r] = append_row(lv2[r], rep_vector(-999, max_cols-dim_2[2]));
	i += dim_2[1];
}

if(dim_3[1] > 0) {
	for(r in 1:dim_3[1]) merged[i+r] = append_row(lv3[r], rep_vector(-999, max_cols-dim_3[2]));
	i += dim_3[1];
}

if(dim_4[1] > 0) {
	for(r in 1:dim_4[1]) merged[i+r] = append_row(lv4[r], rep_vector(-999, max_cols-dim_4[2]));
}

	return(merged);

}






	int[,] package_int(
		int C,
		int GM,
		int Q,
		int shards,
		int[] counts,

		int[] size_counts_idx_lv_1_MPI,
		int[] size_counts_G_lv_1_MPI,
		int[] size_counts_G_lv_1_MPI_non_redundant,
		int[] size_counts_S_lv_1_MPI,
		int[] size_G_linear_MPI,
		int[] size_y_linear_S_1_MPI,
		int[] size_y_linear_1_MPI,
		int[,] counts_idx_lv_1_MPI,
		int[,] counts_G_lv_1_MPI_non_redundant_reps,

		int[,]	y_linear_1_MPI

	){
		int threshold = -999;
		int dim_indices = 10;

		int max_col = dim_indices + max(size_counts_idx_lv_1_MPI)+ max(size_counts_G_lv_1_MPI_non_redundant) + max(size_y_linear_1_MPI) ;

		int int_pack[shards,max_col];

		for(i in 1:shards){

			int real_col = dim_indices + size_counts_idx_lv_1_MPI[i]+  size_counts_G_lv_1_MPI_non_redundant[i] + size_y_linear_1_MPI[i] ;

			int_pack[i] =
				concatenate_int_array({
					// indexes
					{C},
					{GM},
					{Q},
					{size_counts_idx_lv_1_MPI[i]},
	 				{size_counts_G_lv_1_MPI[i]},
	 				{size_counts_G_lv_1_MPI_non_redundant[i]},
		 			{size_counts_S_lv_1_MPI[i]},
		 			{size_G_linear_MPI[i]},
		 			{size_y_linear_S_1_MPI[i]},
		 			{size_y_linear_1_MPI[i]},

					// data
					counts[counts_idx_lv_1_MPI[i, 1:size_counts_idx_lv_1_MPI[i]]],

					// Indexed for optimised MPI
					counts_G_lv_1_MPI_non_redundant_reps[i, 1:size_counts_G_lv_1_MPI_non_redundant[i]],

					// counts mix
					y_linear_1_MPI[i, 1:size_y_linear_1_MPI[i]],


					// Buffer
					rep_int(threshold, max_col-real_col)
				});

			//print("|| int ", num_elements(int_pack[i]));
		}

		return int_pack;
	}

	vector[] package_real(
		int C,
		int Q,
		int shards,

		int[] size_counts_G_lv_1_MPI,
		int[] size_counts_G_lv_1_MPI_non_redundant,
		int[] size_counts_S_lv_1_MPI,
		int[] size_G_linear_MPI,
		int[] size_y_linear_S_1_MPI,

		int[,] counts_G_lv_1_MPI,
		int[,] counts_G_lv_1_MPI_non_redundant,
		int[,] counts_S_lv_1_MPI,

		int[,] G_linear_MPI,
		int[,] y_linear_S_1_MPI,

		vector lambda_log,
		vector exposure_rate,
		vector sigma_inv_log,
		vector[] prop_1
	){
		real threshold = -999;

		int max_col =
			max(size_counts_G_lv_1_MPI_non_redundant) +
			max(size_counts_S_lv_1_MPI) +
			max(size_counts_G_lv_1_MPI_non_redundant) +
			(Q*C) +
			max(size_G_linear_MPI) +
			max(size_G_linear_MPI) +
			max(size_y_linear_S_1_MPI);

		vector[max_col] real_pack[shards];

		for(i in 1:shards){

			int real_col =
				size_counts_G_lv_1_MPI_non_redundant[i] +
				size_counts_S_lv_1_MPI[i] +
				size_counts_G_lv_1_MPI_non_redundant[i] +
				(Q*C) +
				size_G_linear_MPI[i] +
				size_G_linear_MPI[i] +
				size_y_linear_S_1_MPI[i];

			real_pack[i] =
				concatenate_vector_array({

					// The gene sample specific redundant estimates of the counts for the reference
					lambda_log[counts_G_lv_1_MPI_non_redundant[i, 1:size_counts_G_lv_1_MPI_non_redundant[i]]],

					// The exposure rate of all samples
					exposure_rate[counts_S_lv_1_MPI[i,1:size_counts_S_lv_1_MPI[i]]],

					// The gene sample specific redundant overdispersions of the counts for the reference
					sigma_inv_log[counts_G_lv_1_MPI_non_redundant[i, 1:size_counts_G_lv_1_MPI_non_redundant[i]]],

					// Proportion vector
					to_vector(vector_array_to_matrix(prop_1)),

					// The estimated of ref to de convolved together
					lambda_log[G_linear_MPI[i, 1:size_G_linear_MPI[i]]],

					// The estimated of ref to de convolved together
					sigma_inv_log[G_linear_MPI[i, 1:size_G_linear_MPI[i]]],

					// The exposure of the mix query samples
					exposure_rate[y_linear_S_1_MPI[i, 1:size_y_linear_S_1_MPI[i]]],

					// Buffer
					rep_vector(threshold, max_col-real_col)
				});

			//print("|| real ", num_elements(real_pack[i]));

		}

		return real_pack;
	}



  vector[] sum_NB_MPI(matrix lambda_mat, matrix sigma_mat, matrix prop_mat){

		// Matrix operation for sum
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

		vector[rows(prop_mat) * cols(sigma_mat)] sum_obj[2];
		sum_obj[1] = to_vector(lambda_sum);
		sum_obj[2] = to_vector(sigma_sum);

		// The vectorisation is G1-Q1, G1-Q2, G1-Q3 etc..
		return(sum_obj);
	}

	vector lp_reduce( vector global_parameters , vector local_parameters , real[] real_data , int[] int_data ) {

		real lp;
		real threshold = -999;

		// int unpacking

		int C = int_data[1];
		int GM = int_data[2];
		int Q = int_data[3];
		int size_counts_idx_lv_1_MPI= int_data[4];
	 	int size_counts_G_lv_1_MPI= int_data[5];
	 	int size_counts_G_lv_1_MPI_non_redundant = int_data[6];
		int size_counts_S_lv_1_MPI= int_data[7];
		int size_G1_linear_MPI= int_data[8];
		int size_y_linear_S_1_MPI= int_data[9];
		int size_y_linear_1_MPI= int_data[10];
		int dim_indices = 10;

		int ref_counts[size_counts_idx_lv_1_MPI] = int_data[dim_indices+1 : dim_indices+size_counts_idx_lv_1_MPI];
		int counts_G_lv_1_MPI_non_redundant_reps[size_counts_G_lv_1_MPI_non_redundant] = int_data[dim_indices+size_counts_idx_lv_1_MPI + 1 : dim_indices + size_counts_idx_lv_1_MPI + size_counts_G_lv_1_MPI_non_redundant];
		int mix_counts[size_y_linear_1_MPI] = int_data[dim_indices+size_counts_idx_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant +1 : dim_indices+size_counts_idx_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+ size_y_linear_1_MPI];

		// real unpacking
		vector[size_counts_G_lv_1_MPI_non_redundant] ref_lambda_log_for_counts = local_parameters[1 : size_counts_G_lv_1_MPI_non_redundant];
		vector[size_counts_S_lv_1_MPI] ref_exposure_rate = local_parameters[size_counts_G_lv_1_MPI_non_redundant+1 : size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI];
		vector[size_counts_G_lv_1_MPI_non_redundant] ref_sigma_inv_log_for_counts =  local_parameters[  size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+1 :  size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant];
		vector[C*Q] prop_1 = local_parameters[size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+1 : size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)];

		vector[size_G1_linear_MPI] ref_lambda_log = local_parameters[size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+1 :size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+size_G1_linear_MPI];

		vector[size_G1_linear_MPI] ref_sigma_inv_log = local_parameters[size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+size_G1_linear_MPI+1 : size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+size_G1_linear_MPI+size_G1_linear_MPI];

		vector[size_y_linear_S_1_MPI] mix_exposure_rate = local_parameters[size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+size_G1_linear_MPI+size_G1_linear_MPI+1 : size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+size_G1_linear_MPI+size_G1_linear_MPI+size_y_linear_S_1_MPI];

		// Ref variables
		vector[size_counts_G_lv_1_MPI] ref_lambda_log_redundant_normalised = rep_vector_by_array(ref_lambda_log_for_counts, counts_G_lv_1_MPI_non_redundant_reps) + ref_exposure_rate;
		vector[size_counts_G_lv_1_MPI] ref_sigma_inv_log_redundant =  rep_vector_by_array(ref_sigma_inv_log_for_counts, counts_G_lv_1_MPI_non_redundant_reps);


		// Mix

		// global parameters
		real sigma_slope = global_parameters[2];
		real sigma_intercept_dec = global_parameters[1];

		// Deconvoluted means
		vector[size_G1_linear_MPI/C*Q] lambda_log_deconvoluted_1;
    vector[size_G1_linear_MPI/C*Q] sigma_deconvoluted_1;

		// Calculate convoluted
    vector[size_G1_linear_MPI/C * Q] sumNB[2] = sum_NB_MPI(
      to_matrix(exp(ref_lambda_log), C, size_G1_linear_MPI/C),
      to_matrix(1.0 ./ exp( ref_sigma_inv_log ), C, size_G1_linear_MPI/C),
      to_matrix(prop_1, Q, C)
    );




    lambda_log_deconvoluted_1 = log(sumNB[1]);
    sigma_deconvoluted_1 = sumNB[2];


// print(mix_counts);
// print(lambda_log_deconvoluted_1);
// print(mix_exposure_rate);
// print(sigma_deconvoluted_1);


		// Reference
		lp = neg_binomial_2_log_lpmf( ref_counts |	ref_lambda_log_redundant_normalised, 1.0 ./ exp( ref_sigma_inv_log_redundant )		);
		//if(sum(to_matrix(prop_1, S, C)[1,]) !=1.0 ) stop("Mpi fucked");

		// deconvolution
		lp += neg_binomial_2_log_lpmf(mix_counts |
			lambda_log_deconvoluted_1 + mix_exposure_rate,
			sigma_deconvoluted_1 //(1.0 ./ exp( ( sigma_slope * lambda_log_deconvoluted_1 + sigma_intercept_dec ) ))
		);


	 return [lp]';

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

	real dirichlet_regression_lpdf(vector p, row_vector X, matrix alpha, real phi){

		// Build sum to zero variable
		int c = cols(alpha);
		int r = rows(alpha);
		matrix[r, c]  alpha_ = alpha;
		alpha_[1,c] = -sum(alpha_[1, 1:(c-1)]);


		// Calculate log prob
		return (dirichlet_lpdf(p | softmax( to_vector(X * alpha_ ))  * exp(phi) + 1 ));
	}

}
data {
	// shards
	int<lower=1> shards;
	int lv;
	// Reference matrix inference
	int<lower=0> G;
	int<lower=0> GM;
	int<lower=0> S;
	int CL; // counts linear size

	// reference counts
 	int<lower=0> counts_linear[CL] ;

	int G_to_counts_linear[CL] ;
	int S_linear[CL] ;

	// Reference counts per level
	int<lower=0> CL_NA;
	int<lower=0> counts_idx_lv_NA[CL_NA];
	int<lower=0> CL_lv;
	int<lower=0> counts_idx_lv[CL_lv];

	// Priors
	real<upper=0> sigma_slope;
	real<lower=0> sigma_sigma;
	real sigma_intercept;

	// Cell types
 	int<lower=0> Q;
  int<lower=1> n_nodes;
  int<lower=1> ct_in_nodes[n_nodes];
	int<lower=1> n_levels;
  int<lower=1> ct_in_levels[n_levels];

  // Deconvolution
  int<lower=0> G_lv;
  int G_lv_linear[G_lv];

  // Observed counts
  int<lower=0> Y_lv;
	int y_linear_lv[Y_lv];
	int y_linear_S_lv[Y_lv];

	// MPI
	int<lower=0> size_counts_idx_lv_MPI[max(shards, 1)];
	int counts_idx_lv_MPI[shards, max(size_counts_idx_lv_MPI)];
	int size_counts_G_lv_MPI[max(shards, 1)];
	int size_counts_S_lv_MPI[max(shards, 1)];
	int size_y_linear_S_MPI[max(shards, 1)];
	int size_y_linear_MPI[max(shards, 1)];
	int y_linear_MPI[shards,max(size_y_linear_MPI)];
	int y_linear_S_MPI[shards,max(size_y_linear_S_MPI)];
	int counts_G_lv_MPI[shards,max(size_counts_G_lv_MPI)];
	int counts_S_lv_MPI[shards,max(size_counts_S_lv_MPI)];
	int size_G_linear_MPI[max(shards, 1)];
	int G_linear_MPI[shards,max(size_G_linear_MPI)];
	int size_counts_G_lv_MPI_non_redundant[max(shards, 1)];
	int counts_G_lv_MPI_non_redundant[shards, max(size_counts_G_lv_MPI_non_redundant)];
	int counts_G_lv_MPI_non_redundant_reps[shards, max(size_counts_G_lv_MPI_non_redundant)];

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

	// Non-centered param
	real lambda_mu_prior[2];
	real lambda_sigma_prior[2];
	real lambda_skew_prior[2];
	real sigma_intercept_prior[2];

  // Proportions priors
  // lv1
  vector[ct_in_nodes[1]]  prop_1_prior[Q * (lv > 1)]; // Root

  // lv2
  vector[ct_in_nodes[2]]  prop_a_prior[Q * (lv > 2)]; // Immune cells

  // lv3
  vector[ct_in_nodes[3]]  prop_b_prior[Q * (lv > 3)]; // b cells
  vector[ct_in_nodes[4]]  prop_c_prior[Q * (lv > 3)]; // granulocyte
  vector[ct_in_nodes[5]]  prop_d_prior[Q * (lv > 3)]; // mono_derived
  vector[ct_in_nodes[6]]  prop_e_prior[Q * (lv > 3)]; // t_cell

	// exposure posterior previous fit
	// vector[2] exposure_posterior[Q * (lv > 1)];
	// real exposure_rate_multiplier;

	// Dirichlet regression
	int A; // factors of interest
	matrix[Q,A] X;
	int do_regression;

}
transformed data{

	real real_data[shards, 0];
	real real_data2[shards, 0];

	int dim_indices = 10;

	int pack_cols = dim_indices + max(size_counts_idx_lv_MPI) + max(size_counts_G_lv_MPI_non_redundant) + max(size_y_linear_MPI);

	int pack[shards, pack_cols];

	int int_package[shards, max({pack_cols })];

	pack = package_int(
		ct_in_levels[lv],
		G_lv/ct_in_levels[lv],
		Q,
		shards,
		counts_linear,
		size_counts_idx_lv_MPI,
		size_counts_G_lv_MPI,
		size_counts_G_lv_MPI_non_redundant,
		size_counts_S_lv_MPI,
		size_G_linear_MPI,
		size_y_linear_S_MPI,
		size_y_linear_MPI,
		counts_idx_lv_MPI,
		counts_G_lv_MPI_non_redundant_reps,
		y_linear_MPI
	);

	// Here I am building the whole int package to not have to calculate it every time
	int_package	= pack	;

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
  // lv1
  simplex[ct_in_nodes[1]]  prop_1[Q * (lv >= 1)]; // Root

  // lv2
  simplex[ct_in_nodes[2]]  prop_a[Q * (lv >= 2)]; // Immune cells

  // lv3
  simplex[ct_in_nodes[3]]  prop_b[Q * (lv >= 3)]; // b cells
  simplex[ct_in_nodes[4]]  prop_c[Q * (lv >= 3)]; // granulocyte
  simplex[ct_in_nodes[5]]  prop_d[Q * (lv >= 3)]; // mono_derived
  simplex[ct_in_nodes[6]]  prop_e[Q * (lv >= 3)]; // t_cell

	// lv4
  simplex[ct_in_nodes[7]]  prop_f[Q * (lv >= 4)]; // dendritic myeloid
  simplex[ct_in_nodes[8]]  prop_g[Q * (lv >= 4)]; // macrophage
  simplex[ct_in_nodes[9]]  prop_h[Q * (lv >= 4)]; // CD4
  simplex[ct_in_nodes[10]] prop_i[Q * (lv >= 4)]; // CD8

  // Dirichlet regression
  // lv1
  matrix[A * (lv == 1) * do_regression,ct_in_nodes[1]]  alpha_1; // Root

  // lv2
  matrix[A * (lv == 2) * do_regression,ct_in_nodes[2]]  alpha_a; // Immune cells

  // lv3
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[3]]  alpha_b; // b cells
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[4]]  alpha_c; // granulocyte
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[5]]  alpha_d; // mono_derived
  matrix[A * (lv == 3) * do_regression,ct_in_nodes[6]]  alpha_e; // t_cell

	// lv4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[7]]  alpha_f; // dendritic myeloid
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[8]]  alpha_g; // macrophage
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[9]]  alpha_h; // CD4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[10]] alpha_i; // CD8

	real phi[10];

}
transformed parameters{

	vector[ct_in_levels[2]] prop_2[Q * (lv >= 2)];
	vector[ct_in_levels[3]] prop_3[Q * (lv >= 3)];
	vector[ct_in_levels[4]] prop_4[Q * (lv >= 4)];

	// proportion of level 2
	if(lv >= 2)
	prop_2 =
		append_vector_array(
			prop_1[,singles_lv2],
			multiply_by_column(prop_a, prop_1[,parents_lv2[1]])
		);

	// proportion of level 3
	if(lv >= 3)
	prop_3 =
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

	// proportion of level 4
	if(lv >= 4)
	prop_4 =
		append_vector_array(
			prop_3[,singles_lv4],
			append_vector_array(
				multiply_by_column(prop_f, prop_3[,parents_lv4[1]]),
				append_vector_array(
					multiply_by_column(prop_g, prop_3[,parents_lv4[2]]),
					append_vector_array(
						multiply_by_column(prop_h, prop_3[,parents_lv4[3]]),
						multiply_by_column(prop_i, prop_3[,parents_lv4[4]])
					)
				)
			)
		);

}
model {

	vector[Y_lv] lambda_log_deconvoluted;
	vector[Y_lv] sigma_deconvoluted;

 	vector[ct_in_levels[lv]] prop_lv[Q] = which(lv, prop_1, prop_2, prop_3, prop_4);

	vector[Y_lv] sumNB[2];

	vector[
		max(size_counts_G_lv_MPI_non_redundant) +
		max(size_counts_S_lv_MPI) +
		max(size_counts_G_lv_MPI_non_redundant) +
		(ct_in_levels[lv] * Q) +
		max(size_G_linear_MPI) +
		max(size_G_linear_MPI) +
		max(size_y_linear_S_MPI)
	] pack_r[shards];

	real lp = 0;
	real lp_MPI = 0;

  // Overall properties of the data
  lambda_mu ~ normal(lambda_mu_prior[1],lambda_mu_prior[2]);
	lambda_sigma ~ normal(lambda_sigma_prior[1],lambda_sigma_prior[2]);
	lambda_skew ~ normal(lambda_skew_prior[1],lambda_skew_prior[2]);

	// Exposure
	exposure_rate ~ normal(0,1);
	sum(exposure_rate) ~ normal(0, 0.001 * S);

	// Means overdispersion reference
	lambda_log ~ skew_normal(lambda_mu, exp(lambda_sigma), lambda_skew);
	sigma_inv_log ~ normal(sigma_slope * lambda_log + sigma_intercept, sigma_sigma);

	// Deconvolution
	sigma_intercept_dec ~ student_t(3, 0, 2);

	// Level NA - Mix house keeing /////////////////////

	target += sum(map_rect(
		lp_reduce_simple ,
		[sigma_intercept, sigma_slope]', // global parameters
		get_mu_sigma_vector_MPI(
			lambda_log[G_to_counts_linear[counts_idx_lv_NA]] + exposure_rate[S_linear[counts_idx_lv_NA]],
			sigma_inv_log[G_to_counts_linear[counts_idx_lv_NA]],
			shards
		),
		real_data,
		get_int_MPI( counts_linear[counts_idx_lv_NA], shards)
	));


	// lv 1
  if(lv == 1 && do_regression) {

  	for(q in 1:Q) prop_1[q] ~ dirichlet_regression( X[q], alpha_1, phi[1] );
  	to_vector( alpha_1 ) ~ normal(0,1);

  }
	if(lv == 1 && !do_regression) for(q in 1:Q) target += dirichlet_lpdf(prop_1[q] | rep_vector(num_elements(prop_1[1]), num_elements(prop_1[1])));
	if(lv > 1)  for(q in 1:Q) target += dirichlet_lpdf(prop_1[q] | prop_1_prior[q]);

	// lv 2
  if(lv == 2 && do_regression) {

  	for(q in 1:Q) prop_a[q] ~ dirichlet_regression( X[q], alpha_a, phi[2] );
  	to_vector( alpha_a ) ~ normal(0,1);

  }
	if(lv == 2 && !do_regression) for(q in 1:Q) target += dirichlet_lpdf(prop_a[q] | rep_vector(num_elements(prop_a[1]), num_elements(prop_a[1])));
	if(lv > 2)  for(q in 1:Q) target += dirichlet_lpdf(prop_a[q] | prop_a_prior[q]);

	// lv 3
  if(lv == 3 && do_regression){

  	for(q in 1:Q) prop_b[q] ~ dirichlet_regression( X[q], alpha_b, phi[3] );
  	to_vector( alpha_b ) ~ normal(0,1);

  	for(q in 1:Q) prop_c[q] ~ dirichlet_regression( X[q], alpha_c, phi[4] );
  	to_vector( alpha_c ) ~ normal(0,1);

  	for(q in 1:Q) prop_d[q] ~ dirichlet_regression( X[q], alpha_d, phi[5] );
  	to_vector( alpha_d ) ~ normal(0,1);

  	for(q in 1:Q) prop_e[q] ~ dirichlet_regression( X[q], alpha_e, phi[6] );
  	to_vector( alpha_e ) ~ normal(0,1);

  }
  if(lv == 3 && !do_regression) for(q in 1:Q){
  	 target += dirichlet_lpdf(prop_b[q] | rep_vector(num_elements(prop_b[1]), num_elements(prop_b[1])));
		 target += dirichlet_lpdf(prop_c[q] | rep_vector(num_elements(prop_c[1]), num_elements(prop_c[1])));
		 target += dirichlet_lpdf(prop_d[q] | rep_vector(num_elements(prop_d[1]), num_elements(prop_d[1])));
		 target += dirichlet_lpdf(prop_e[q] | rep_vector(num_elements(prop_e[1]), num_elements(prop_e[1])));
  }

	// lv 4
	if(lv == 4) for(q in 1:Q) {
		target += dirichlet_lpdf(prop_f[q] | rep_vector(num_elements(prop_f[1]), num_elements(prop_f[1])));
		target += dirichlet_lpdf(prop_g[q] | rep_vector(num_elements(prop_g[1]), num_elements(prop_g[1])));
		target += dirichlet_lpdf(prop_h[q] | rep_vector(num_elements(prop_h[1]), num_elements(prop_h[1])));
		target += dirichlet_lpdf(prop_i[q] | rep_vector(num_elements(prop_i[1]), num_elements(prop_i[1])));
	}


	pack_r = package_real(
		ct_in_levels[lv],
		Q,
		shards,
		size_counts_G_lv_MPI,
		size_counts_G_lv_MPI_non_redundant,
		size_counts_S_lv_MPI,
		size_G_linear_MPI,
		size_y_linear_S_MPI,
		counts_G_lv_MPI,
		counts_G_lv_MPI_non_redundant,
		counts_S_lv_MPI,
		G_linear_MPI,
		y_linear_S_MPI,
		lambda_log,
		exposure_rate,
		sigma_inv_log,
		prop_lv
	);

	target +=  sum(map_rect(
		lp_reduce ,
		[sigma_intercept_dec, sigma_slope]' ,
		pack_r,
		real_data2,
		int_package
	));

	// Dirichlet regression
	phi ~ normal(0,1);

}
