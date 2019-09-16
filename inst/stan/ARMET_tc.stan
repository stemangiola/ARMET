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
		int[] size_G1_linear_MPI,
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
		 			{size_G1_linear_MPI[i]},
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
		int[] size_G1_linear_MPI,
		int[] size_y_linear_S_1_MPI,

		int[,] counts_G_lv_1_MPI,
		int[,] counts_G_lv_1_MPI_non_redundant,
		int[,] counts_S_lv_1_MPI,

		int[,] G1_linear_MPI,
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
			max(size_G1_linear_MPI) +
			max(size_y_linear_S_1_MPI);

		vector[max_col] real_pack[shards];

		for(i in 1:shards){

			int real_col =
				size_counts_G_lv_1_MPI_non_redundant[i] +
				size_counts_S_lv_1_MPI[i] +
				size_counts_G_lv_1_MPI_non_redundant[i] +
				(Q*C) +
				size_G1_linear_MPI[i] +
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
					lambda_log[G1_linear_MPI[i, 1:size_G1_linear_MPI[i]]],

					// The exposure of the mix query samples
					exposure_rate[y_linear_S_1_MPI[i, 1:size_y_linear_S_1_MPI[i]]],

					// Buffer
					rep_vector(threshold, max_col-real_col)
				});

			//print("|| real ", num_elements(real_pack[i]));

		}

		return real_pack;
	}


	int[,] append_int_MPI_arrays(int[,] lv1, int[,] lv2, int[,] lv3, int[,] lv4){
		int dim_1[2] = dims(lv1);
		int dim_2[2] = dims(lv2);
		int dim_3[2] = dims(lv3);
		int dim_4[2] = dims(lv4);

		int max_cols = max({dim_1[2], dim_2[2],dim_3[2],dim_4[2]});
		int tot_rows = dim_1[1] + dim_2[1] + dim_3[1] + dim_4[1];


		int merged[tot_rows, max_cols];

		int i = 0;


		for(r in 1:dim_1[1]) merged[i+r] = append_array(lv1[r], rep_int(-999, max_cols-dim_1[2]));

		i += dim_1[1];
		for(r in 1:dim_2[1]) merged[i+r] = append_array(lv2[r], rep_int(-999, max_cols-dim_2[2]));

		i += dim_2[1];
		for(r in 1:dim_3[1]) merged[i+r] = append_array(lv3[r], rep_int(-999, max_cols-dim_3[2]));

		i += dim_3[1];
		for(r in 1:dim_4[1]) merged[i+r] = append_array(lv4[r], rep_int(-999, max_cols-dim_4[2]));


		return(merged);

	}

	vector[] append_real_MPI_arrays(vector[] lv1, vector[] lv2, vector[] lv3, vector[] lv4){
		int dim_1[2] = {num_elements(lv1[,1]), num_elements(lv1[1])};
		int dim_2[2] = {num_elements(lv2[,1]), num_elements(lv2[1])};
		int dim_3[2] = {num_elements(lv3[,1]), num_elements(lv3[1])};
		int dim_4[2] = {num_elements(lv4[,1]), num_elements(lv4[1])};

		int max_cols = max({dim_1[2], dim_2[2],dim_3[2],dim_4[2]});
		int tot_rows = dim_1[1] + dim_2[1] + dim_3[1] + dim_4[1];

		vector[max_cols] merged[tot_rows];
		int i = 0;


		for(r in 1:dim_1[1]) merged[i+r] = append_row(lv1[r], rep_vector(-999, max_cols-dim_1[2]));

		i += dim_1[1];
		for(r in 1:dim_2[1]) merged[i+r] = append_row(lv2[r], rep_vector(-999, max_cols-dim_2[2]));

		i += dim_2[1];
		for(r in 1:dim_3[1]) merged[i+r] = append_row(lv3[r], rep_vector(-999, max_cols-dim_3[2]));

		i += dim_3[1];
		for(r in 1:dim_4[1]) merged[i+r] = append_row(lv4[r], rep_vector(-999, max_cols-dim_4[2]));


		return(merged);

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
		vector[size_counts_G_lv_1_MPI_non_redundant] ref_sigma_inv_log =  local_parameters[  size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+1 :  size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant];
		vector[C*Q] prop_1 = local_parameters[size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+1 : size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)];
		vector[size_G1_linear_MPI] ref_lambda_log = local_parameters[size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+1 :size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+size_G1_linear_MPI];
		vector[size_y_linear_S_1_MPI] mix_exposure_rate = local_parameters[size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+size_G1_linear_MPI+1 : size_counts_G_lv_1_MPI_non_redundant+size_counts_S_lv_1_MPI+size_counts_G_lv_1_MPI_non_redundant+(C*Q)+size_G1_linear_MPI+size_y_linear_S_1_MPI];

		// Ref variables
		vector[size_counts_G_lv_1_MPI] ref_lambda_log_redundant_normalised = rep_vector_by_array(ref_lambda_log_for_counts, counts_G_lv_1_MPI_non_redundant_reps) + ref_exposure_rate;
		vector[size_counts_G_lv_1_MPI] ref_sigma_inv_log_redundant =  rep_vector_by_array(ref_sigma_inv_log, counts_G_lv_1_MPI_non_redundant_reps);

		//print(size_ref_counts , "..", size_ref_counts, "..",(C*Q), "..",(GM*C));

		// Mix

		// global parameters
		real sigma_slope = global_parameters[2];
		real sigma_intercept_dec = global_parameters[1];

		// Deconvoluted means
		vector[size_G1_linear_MPI/C*Q] lambda_log_deconvoluted_1;

		// Reference
		lp = neg_binomial_2_log_lpmf( ref_counts |	ref_lambda_log_redundant_normalised, 1.0 ./ exp( ref_sigma_inv_log_redundant )		);

		//print(to_matrix(prop_1, Q, C) );
		//print(size_G1_linear_MPI/C);

		// Calculate convoluted
		lambda_log_deconvoluted_1 =
			log(
				to_vector(
					to_matrix(prop_1, Q, C) *
					exp(to_matrix(ref_lambda_log, C, size_G1_linear_MPI/C)) // [Q,G] dimensions
				)
			);

		//if(sum(to_matrix(prop_1, S, C)[1,]) !=1.0 ) stop("Mpi fucked");

		// deconvolution
		lp += neg_binomial_2_log_lpmf(mix_counts |
			lambda_log_deconvoluted_1 + mix_exposure_rate,
			1.0 ./ exp( ( sigma_slope * lambda_log_deconvoluted_1 + sigma_intercept_dec ) )
		);


	 return [lp]';

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

	// Reference counts per level
	int<lower=0> CL_NA;
	int<lower=0> counts_idx_lv_NA[CL_NA];
	int<lower=0> CL_1;
	int<lower=0> counts_idx_lv_1[CL_1];
	int<lower=0> CL_2;
	int<lower=0> counts_idx_lv_2[CL_2];
	int<lower=0> CL_3;
	int<lower=0> counts_idx_lv_3[CL_3];
	int<lower=0> CL_4;
	int<lower=0> counts_idx_lv_4[CL_4];



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
  int<lower=0> G1;
  int G1_linear[G1];
  int<lower=0> G2;
  int G2_linear[G2];
  int<lower=0> G3;
  int G3_linear[G3];
  int<lower=0> G4;
  int G4_linear[G4];

  // Observed counts
  int<lower=0> Y_1;
	int y_linear_1[Y_1];
	int y_linear_S_1[Y_1];
  int<lower=0> Y_2;
	int y_linear_2[Y_2];
	int y_linear_S_2[Y_2];
	int<lower=0> Y_3;
	int y_linear_3[Y_3];
	int y_linear_S_3[Y_3];
	int<lower=0> Y_4;
	int y_linear_4[Y_4];
	int y_linear_S_4[Y_4];

	// MPI
	int<lower=1> shards_in_levels[4];

	// MPI lv1
	int<lower=0> size_counts_idx_lv_1_MPI[shards_in_levels[1]];
	int counts_idx_lv_1_MPI[shards_in_levels[1], max(size_counts_idx_lv_1_MPI)];
	int size_counts_G_lv_1_MPI[shards_in_levels[1]];
	int size_counts_S_lv_1_MPI[shards_in_levels[1]];
	int size_y_linear_S_1_MPI[shards_in_levels[1]];
	int size_y_linear_1_MPI[shards_in_levels[1]];
	int y_linear_1_MPI[shards_in_levels[1],max(size_y_linear_1_MPI)];
	int y_linear_S_1_MPI[shards_in_levels[1],max(size_y_linear_1_MPI)];
	int counts_G_lv_1_MPI[shards_in_levels[1],max(size_counts_G_lv_1_MPI)];
	int counts_S_lv_1_MPI[shards_in_levels[1],max(size_counts_S_lv_1_MPI)];
	int size_G1_linear_MPI[shards_in_levels[1]];
	int G1_linear_MPI[shards_in_levels[1],max(size_G1_linear_MPI)];
	int size_counts_G_lv_1_MPI_non_redundant[shards_in_levels[1]];
	int counts_G_lv_1_MPI_non_redundant[shards_in_levels[1], max(size_counts_G_lv_1_MPI_non_redundant)];
	int counts_G_lv_1_MPI_non_redundant_reps[shards_in_levels[1], max(size_counts_G_lv_1_MPI_non_redundant)];

	// MPI lv2
	int<lower=0> size_counts_idx_lv_2_MPI[shards_in_levels[2]];
	int counts_idx_lv_2_MPI[shards_in_levels[2], max(size_counts_idx_lv_2_MPI)];
	int size_counts_G_lv_2_MPI[shards_in_levels[2]];
	int size_counts_S_lv_2_MPI[shards_in_levels[2]];
	int size_y_linear_S_2_MPI[shards_in_levels[2]];
	int size_y_linear_2_MPI[shards_in_levels[2]];
	int y_linear_2_MPI[shards_in_levels[2],max(size_y_linear_2_MPI)];
	int y_linear_S_2_MPI[shards_in_levels[2],max(size_y_linear_2_MPI)];
	int counts_G_lv_2_MPI[shards_in_levels[2],max(size_counts_G_lv_2_MPI)];
	int counts_S_lv_2_MPI[shards_in_levels[2],max(size_counts_S_lv_2_MPI)];
	int size_G2_linear_MPI[shards_in_levels[2]];
	int G2_linear_MPI[shards_in_levels[2],max(size_G2_linear_MPI)];
	int size_counts_G_lv_2_MPI_non_redundant[shards_in_levels[2]];
	int counts_G_lv_2_MPI_non_redundant[shards_in_levels[2], max(size_counts_G_lv_2_MPI_non_redundant)];
	int counts_G_lv_2_MPI_non_redundant_reps[shards_in_levels[2], max(size_counts_G_lv_2_MPI_non_redundant)];

	// MPI lv3
	int<lower=0> size_counts_idx_lv_3_MPI[shards_in_levels[3]];
	int counts_idx_lv_3_MPI[shards_in_levels[3], max(size_counts_idx_lv_3_MPI)];
	int size_counts_G_lv_3_MPI[shards_in_levels[3]];
	int size_counts_S_lv_3_MPI[shards_in_levels[3]];
	int size_y_linear_S_3_MPI[shards_in_levels[3]];
	int size_y_linear_3_MPI[shards_in_levels[3]];
	int y_linear_3_MPI[shards_in_levels[3],max(size_y_linear_3_MPI)];
	int y_linear_S_3_MPI[shards_in_levels[3],max(size_y_linear_3_MPI)];
	int counts_G_lv_3_MPI[shards_in_levels[3],max(size_counts_G_lv_3_MPI)];
	int counts_S_lv_3_MPI[shards_in_levels[3],max(size_counts_S_lv_3_MPI)];
	int size_G3_linear_MPI[shards_in_levels[3]];
	int G3_linear_MPI[shards_in_levels[3],max(size_G3_linear_MPI)];
	int size_counts_G_lv_3_MPI_non_redundant[shards_in_levels[3]];
	int counts_G_lv_3_MPI_non_redundant[shards_in_levels[3], max(size_counts_G_lv_3_MPI_non_redundant)];
	int counts_G_lv_3_MPI_non_redundant_reps[shards_in_levels[3], max(size_counts_G_lv_3_MPI_non_redundant)];

	// MPI lv4
	int<lower=0> size_counts_idx_lv_4_MPI[shards_in_levels[4]];
	int counts_idx_lv_4_MPI[shards_in_levels[4], max(size_counts_idx_lv_4_MPI)];
	int size_counts_G_lv_4_MPI[shards_in_levels[4]];
	int size_counts_S_lv_4_MPI[shards_in_levels[4]];
	int size_y_linear_S_4_MPI[shards_in_levels[4]];
	int size_y_linear_4_MPI[shards_in_levels[4]];
	int y_linear_4_MPI[shards_in_levels[4],max(size_y_linear_4_MPI)];
	int y_linear_S_4_MPI[shards_in_levels[4],max(size_y_linear_4_MPI)];
	int counts_G_lv_4_MPI[shards_in_levels[4],max(size_counts_G_lv_4_MPI)];
	int counts_S_lv_4_MPI[shards_in_levels[4],max(size_counts_S_lv_4_MPI)];
	int size_G4_linear_MPI[shards_in_levels[4]];
	int G4_linear_MPI[shards_in_levels[4],max(size_G4_linear_MPI)];
	int size_counts_G_lv_4_MPI_non_redundant[shards_in_levels[4]];
	int counts_G_lv_4_MPI_non_redundant[shards_in_levels[4], max(size_counts_G_lv_4_MPI_non_redundant)];
	int counts_G_lv_4_MPI_non_redundant_reps[shards_in_levels[4], max(size_counts_G_lv_4_MPI_non_redundant)];


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

	// Lv4 tree structure parents singles
	int<lower=1> SLV4;
	int<lower=1> PLV4;
	int parents_lv4[PLV4]; // Level one parents
	int singles_lv4[SLV4]; // Level 1 leafs

	// Non-centered param
	real lambda_mu_prior[2];
	real lambda_sigma_prior[2];
	real lambda_skew_prior[2];
	real sigma_intercept_prior[2];


}
transformed data{

	real real_data[shards, 0];
	real real_data2[sum(shards_in_levels), 0];

	int dim_indices = 10;

	// Here I am building the whole int package to not have to calculate it every time
	int int_package[sum(shards_in_levels), dim_indices +	max({
			max(size_counts_idx_lv_1_MPI) + max(size_counts_G_lv_1_MPI_non_redundant) + max(size_y_linear_1_MPI),
			max(size_counts_idx_lv_2_MPI) + max(size_counts_G_lv_2_MPI_non_redundant) + max(size_y_linear_2_MPI),
			max(size_counts_idx_lv_3_MPI) + max(size_counts_G_lv_3_MPI_non_redundant) + max(size_y_linear_3_MPI),
			max(size_counts_idx_lv_4_MPI) + max(size_counts_G_lv_4_MPI_non_redundant) + max(size_y_linear_4_MPI)
		})] = append_int_MPI_arrays(
			package_int(			ct_in_levels[1],			G1/ct_in_levels[1],			Q,			shards_in_levels[1],			counts_linear,			size_counts_idx_lv_1_MPI,			size_counts_G_lv_1_MPI,	size_counts_G_lv_1_MPI_non_redundant,			size_counts_S_lv_1_MPI,			size_G1_linear_MPI,			size_y_linear_S_1_MPI,			size_y_linear_1_MPI,			counts_idx_lv_1_MPI,	counts_G_lv_1_MPI_non_redundant_reps,		y_linear_1_MPI		),
			package_int(			ct_in_levels[2],			G2/ct_in_levels[2],			Q,			shards_in_levels[2],			counts_linear,			size_counts_idx_lv_2_MPI,			size_counts_G_lv_2_MPI,	size_counts_G_lv_2_MPI_non_redundant,			size_counts_S_lv_2_MPI,			size_G2_linear_MPI,			size_y_linear_S_2_MPI,			size_y_linear_2_MPI,			counts_idx_lv_2_MPI,	counts_G_lv_2_MPI_non_redundant_reps,				y_linear_2_MPI		),
			package_int(			ct_in_levels[3],			G3/ct_in_levels[3],			Q,			shards_in_levels[3],			counts_linear,			size_counts_idx_lv_3_MPI,			size_counts_G_lv_3_MPI,	size_counts_G_lv_3_MPI_non_redundant,			size_counts_S_lv_3_MPI,			size_G3_linear_MPI,			size_y_linear_S_3_MPI,			size_y_linear_3_MPI,			counts_idx_lv_3_MPI,	counts_G_lv_3_MPI_non_redundant_reps,				y_linear_3_MPI		),
			package_int(			ct_in_levels[4],			G4/ct_in_levels[4],			Q,			shards_in_levels[4],			counts_linear,			size_counts_idx_lv_4_MPI,			size_counts_G_lv_4_MPI,	size_counts_G_lv_4_MPI_non_redundant,			size_counts_S_lv_4_MPI,			size_G4_linear_MPI,			size_y_linear_S_4_MPI,			size_y_linear_4_MPI,			counts_idx_lv_4_MPI,	counts_G_lv_4_MPI_non_redundant_reps,				y_linear_4_MPI		)

	);

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
  simplex[ct_in_nodes[1]]  prop_1[Q]; // Root

  // lv2
  simplex[ct_in_nodes[2]]  prop_a[Q]; // Immune cells

  // lv3
  simplex[ct_in_nodes[3]]  prop_b[Q]; // b cells
  simplex[ct_in_nodes[4]]  prop_c[Q]; // granulocyte
  simplex[ct_in_nodes[5]]  prop_d[Q]; // mono_derived
  simplex[ct_in_nodes[6]]  prop_e[Q]; // t_cell

	// lv4
  simplex[ct_in_nodes[7]]  prop_f[Q]; // dendritic myeloid
  simplex[ct_in_nodes[8]]  prop_g[Q]; // macrophage
  simplex[ct_in_nodes[9]]  prop_h[Q]; // CD4
  simplex[ct_in_nodes[10]] prop_i[Q]; // CD8

  // Error between reference and mix, to avoid divergencies
  // vector<lower=0>[GM] error_ref_mix_z;

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

	// proportion of level 4
	vector[ct_in_levels[4]] prop_4[Q] =
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

	// vector[Y_1] lambda_log_deconvoluted_1;
	// vector[Y_2] lambda_log_deconvoluted_2;
	// vector[Y_3] lambda_log_deconvoluted_3;
	// vector[Y_4] lambda_log_deconvoluted_4;

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

	// Reference
	target += neg_binomial_2_log_lpmf( counts_linear[counts_idx_lv_NA] |
		lambda_log[G_to_counts_linear[counts_idx_lv_NA]] + exposure_rate[S_linear[counts_idx_lv_NA]],
		1.0 ./ exp( sigma_inv_log[G_to_counts_linear[counts_idx_lv_NA]] )
	);

	// Level 1 ////////////////////////////////////////


	// // Reference
	// lp += neg_binomial_2_log_lpmf( counts_linear[counts_idx_lv_1] |
	// 	lambda_log[G_to_counts_linear[counts_idx_lv_1]] + exposure_rate[S_linear[counts_idx_lv_1]],
	// 	1.0 ./ exp( sigma_inv_log[G_to_counts_linear[counts_idx_lv_1]] )
	// );
	//
	// // Calculate convoluted
	// lambda_log_deconvoluted_1 =
	// 	log(
	// 		to_vector(
	// 			vector_array_to_matrix(prop_1) *
	// 			exp(to_matrix(lambda_log[G1_linear], ct_in_levels[1], G1/ct_in_levels[1])) // [Q,G] dimensions
	// 		)
	// 	);
	//
	// // deconvolution
	// lp += neg_binomial_2_log_lpmf(y_linear_1 |
	// 	lambda_log_deconvoluted_1 + exposure_rate[y_linear_S_1],
	// 	1.0 ./ exp( ( sigma_slope * lambda_log_deconvoluted_1 + sigma_intercept_dec ) )
	// );


	//Gene-wise properties of the data
	// target += sum(map_rect(
	// 	lp_reduce ,
	// 	[sigma_intercept_dec, sigma_slope]' ,
	// 	package_real(
	// 		ct_in_levels[1],
	// 		Q,
	// 		shards,
	//
	// 		size_counts_G_lv_1_MPI,
	// 		size_counts_S_lv_1_MPI,
	// 		size_G1_linear_MPI,
	// 		size_y_linear_S_1_MPI,
	//
	// 		counts_G_lv_1_MPI,
	// 		counts_S_lv_1_MPI,
	// 		G1_linear_MPI,
	// 		y_linear_S_1_MPI,
	//
	// 		lambda_log,
	// 		exposure_rate,
	// 		sigma_inv_log,
	// 		prop_1
	// 	),
	// 	real_data,
	// 	package_int(
	// 		ct_in_levels[1],
	// 		G1/ct_in_levels[1],
	// 		Q,
	// 		shards,
	//
	// 		counts_linear,
	//
	// 		size_counts_idx_lv_1_MPI,
	// 		size_counts_G_lv_1_MPI,
	// 		size_counts_S_lv_1_MPI,
	// 		size_G1_linear_MPI,
	// 		size_y_linear_S_1_MPI,
	// 		size_y_linear_1_MPI,
	// 		counts_idx_lv_1_MPI,
	// 		y_linear_1_MPI
	//
	// 	)
	//
	// ));


	for(q in 1:Q) target += dirichlet_lpdf(prop_1[q] | rep_vector(num_elements(prop_1[1]), num_elements(prop_1[1])));

// Level 2 ////////////////////////////////////////

	// // Reference
	// lp += neg_binomial_2_log_lpmf(counts_linear[counts_idx_lv_2] |
	// 	lambda_log[G_to_counts_linear[counts_idx_lv_2]] + exposure_rate[S_linear[counts_idx_lv_2]],
	// 	1.0 ./ exp( sigma_inv_log[G_to_counts_linear[counts_idx_lv_2]] )
	// );
	//
	// // Calculate convoluted
	// lambda_log_deconvoluted_2 =
	// 	log(
	// 		to_vector(
	// 			vector_array_to_matrix(prop_2) *
	// 			exp(to_matrix(lambda_log[G2_linear], ct_in_levels[2], G2/ct_in_levels[2])) // [Q,G] dimensions
	// 		)
	// 	);
	//
	// // deconvolution
	// lp += neg_binomial_2_log_lpmf(y_linear_2 |
	// 	lambda_log_deconvoluted_2 + exposure_rate[y_linear_S_2],
	// 	1.0 ./ exp( ( sigma_slope * lambda_log_deconvoluted_2 + sigma_intercept_dec ) )
	// );

	// target += sum(map_rect(
	// 	lp_reduce ,
	// 	[sigma_intercept_dec, sigma_slope]' ,
	// 	package_real(
	// 		ct_in_levels[2],
	// 		Q,
	// 		shards,
	//
	// 		size_counts_G_lv_2_MPI,
	// 		size_counts_S_lv_2_MPI,
	// 		size_G2_linear_MPI,
	// 		size_y_linear_S_2_MPI,
	//
	// 		counts_G_lv_2_MPI,
	// 		counts_S_lv_2_MPI,
	// 		G2_linear_MPI,
	// 		y_linear_S_2_MPI,
	//
	// 		lambda_log,
	// 		exposure_rate,
	// 		sigma_inv_log,
	// 		prop_2
	// 	),
	// 	real_data,
	// 	package_int(
	// 		ct_in_levels[2],
	// 		G2/ct_in_levels[2],
	// 		Q,
	// 		shards,
	//
	// 		counts_linear,
	//
	// 		size_counts_idx_lv_2_MPI,
	// 		size_counts_G_lv_2_MPI,
	// 		size_counts_S_lv_2_MPI,
	// 		size_G2_linear_MPI,
	// 		size_y_linear_S_2_MPI,
	// 		size_y_linear_2_MPI,
	// 		counts_idx_lv_2_MPI,
	// 		y_linear_2_MPI
	//
	// 	)
	//
	// ));

	for(q in 1:Q) target += dirichlet_lpdf(prop_a[q] | rep_vector(num_elements(prop_a[1]), num_elements(prop_a[1])));

// Level 3 ////////////////////////////////////////

	// // Reference
	// lp += neg_binomial_2_log_lpmf(counts_linear[counts_idx_lv_3] |
	// 	lambda_log[G_to_counts_linear[counts_idx_lv_3]] + exposure_rate[S_linear[counts_idx_lv_3]],
	// 	1.0 ./ exp( sigma_inv_log[G_to_counts_linear[counts_idx_lv_3]] )
	// );
	//
	// // Calculate convoluted
	// lambda_log_deconvoluted_3 =
	// 	log(
	// 		to_vector(
	// 			vector_array_to_matrix(prop_3) *
	// 			exp(to_matrix(lambda_log[G3_linear], ct_in_levels[3], G3/ct_in_levels[3])) // [Q,G] dimensions
	// 		)
	// 	);
	//
	// // deconvolution
	// lp += neg_binomial_2_log_lpmf(y_linear_3 |
	// 	lambda_log_deconvoluted_3 + exposure_rate[y_linear_S_3],
	// 	1.0 ./ exp( ( sigma_slope * lambda_log_deconvoluted_3 + sigma_intercept_dec ) )
	// );

	// target += sum(map_rect(
	// 	lp_reduce ,
	// 	[sigma_intercept_dec, sigma_slope]' ,
	// 	package_real(
	// 		ct_in_levels[3],
	// 		Q,
	// 		shards,
	//
	// 		size_counts_G_lv_3_MPI,
	// 		size_counts_S_lv_3_MPI,
	// 		size_G3_linear_MPI,
	// 		size_y_linear_S_3_MPI,
	//
	// 		counts_G_lv_3_MPI,
	// 		counts_S_lv_3_MPI,
	// 		G3_linear_MPI,
	// 		y_linear_S_3_MPI,
	//
	// 		lambda_log,
	// 		exposure_rate,
	// 		sigma_inv_log,
	// 		prop_3
	// 	),
	// 	real_data,
	// 	package_int(
	// 		ct_in_levels[3],
	// 		G3/ct_in_levels[3],
	// 		Q,
	// 		shards,
	//
	// 		counts_linear,
	//
	// 		size_counts_idx_lv_3_MPI,
	// 		size_counts_G_lv_3_MPI,
	// 		size_counts_S_lv_3_MPI,
	// 		size_G3_linear_MPI,
	// 		size_y_linear_S_3_MPI,
	// 		size_y_linear_3_MPI,
	// 		counts_idx_lv_3_MPI,
	// 		y_linear_3_MPI
	//
	// 	)
	//
	// ));

	for(q in 1:Q) target += dirichlet_lpdf(prop_b[q] | rep_vector(num_elements(prop_b[1]), num_elements(prop_b[1])));
	for(q in 1:Q) target += dirichlet_lpdf(prop_c[q] | rep_vector(num_elements(prop_c[1]), num_elements(prop_c[1])));
	for(q in 1:Q) target += dirichlet_lpdf(prop_d[q] | rep_vector(num_elements(prop_d[1]), num_elements(prop_d[1])));
	for(q in 1:Q) target += dirichlet_lpdf(prop_e[q] | rep_vector(num_elements(prop_e[1]), num_elements(prop_e[1])));

// Level 4 ////////////////////////////////////////

	// // Reference
	// lp += neg_binomial_2_log_lpmf( counts_linear[counts_idx_lv_4] |
	// 	lambda_log[G_to_counts_linear[counts_idx_lv_4]] + exposure_rate[S_linear[counts_idx_lv_4]],
	// 	1.0 ./ exp( sigma_inv_log[G_to_counts_linear[counts_idx_lv_4]] )
	// );
	//
	// // Calculate convoluted
	// lambda_log_deconvoluted_4 =
	// 	log(
	// 		to_vector(
	// 			vector_array_to_matrix(prop_4) *
	// 			exp(to_matrix(lambda_log[G4_linear], ct_in_levels[4], G4/ct_in_levels[4])) // [Q,G] dimensions
	// 		)
	// 	);
	//
	// // deconvolution
	// lp += neg_binomial_2_log_lpmf(y_linear_4 |
	// 	lambda_log_deconvoluted_4 + exposure_rate[y_linear_S_4],
	// 	1.0 ./ exp( ( sigma_slope * lambda_log_deconvoluted_4 + sigma_intercept_dec ) )
	// );

	// target += sum(map_rect(
	// 	lp_reduce ,
	// 	[sigma_intercept_dec, sigma_slope]' ,
	// 	package_real(
	// 		ct_in_levels[4],
	// 		Q,
	// 		shards,
	//
	// 		size_counts_G_lv_4_MPI,
	// 		size_counts_S_lv_4_MPI,
	// 		size_G4_linear_MPI,
	// 		size_y_linear_S_4_MPI,
	//
	// 		counts_G_lv_4_MPI,
	// 		counts_S_lv_4_MPI,
	// 		G4_linear_MPI,
	// 		y_linear_S_4_MPI,
	//
	// 		lambda_log,
	// 		exposure_rate,
	// 		sigma_inv_log,
	// 		prop_4
	// 	),
	// 	real_data,
	// 	package_int(
	// 		ct_in_levels[4],
	// 		G4/ct_in_levels[4],
	// 		Q,
	// 		shards,
	//
	// 		counts_linear,
	//
	// 		size_counts_idx_lv_4_MPI,
	// 		size_counts_G_lv_4_MPI,
	// 		size_counts_S_lv_4_MPI,
	// 		size_G4_linear_MPI,
	// 		size_y_linear_S_4_MPI,
	// 		size_y_linear_4_MPI,
	// 		counts_idx_lv_4_MPI,
	// 		y_linear_4_MPI
	//
	// 	)
	//
	// ));

	for(q in 1:Q) target += dirichlet_lpdf(prop_f[q] | rep_vector(num_elements(prop_f[1]), num_elements(prop_f[1])));
	for(q in 1:Q) target += dirichlet_lpdf(prop_g[q] | rep_vector(num_elements(prop_g[1]), num_elements(prop_g[1])));
	for(q in 1:Q) target += dirichlet_lpdf(prop_h[q] | rep_vector(num_elements(prop_h[1]), num_elements(prop_h[1])));
	for(q in 1:Q) target += dirichlet_lpdf(prop_i[q] | rep_vector(num_elements(prop_i[1]), num_elements(prop_i[1])));


	target += sum(map_rect(
		lp_reduce ,
		[sigma_intercept_dec, sigma_slope]' ,
		append_real_MPI_arrays(
					package_real(	ct_in_levels[1],	Q,	shards_in_levels[1],	size_counts_G_lv_1_MPI, size_counts_G_lv_1_MPI_non_redundant,	size_counts_S_lv_1_MPI,	size_G1_linear_MPI,	size_y_linear_S_1_MPI, counts_G_lv_1_MPI, counts_G_lv_1_MPI_non_redundant, counts_S_lv_1_MPI,		G1_linear_MPI,			y_linear_S_1_MPI,			lambda_log,			exposure_rate,			sigma_inv_log,			prop_1		),
					package_real(	ct_in_levels[2],	Q,	shards_in_levels[2],	size_counts_G_lv_2_MPI,	size_counts_G_lv_2_MPI_non_redundant,	size_counts_S_lv_2_MPI,	size_G2_linear_MPI,	size_y_linear_S_2_MPI, counts_G_lv_2_MPI, counts_G_lv_2_MPI_non_redundant, counts_S_lv_2_MPI,		G2_linear_MPI,			y_linear_S_2_MPI,			lambda_log,			exposure_rate,			sigma_inv_log,			prop_2		),
					package_real(	ct_in_levels[3],	Q,	shards_in_levels[3],	size_counts_G_lv_3_MPI,	size_counts_G_lv_3_MPI_non_redundant,	size_counts_S_lv_3_MPI,	size_G3_linear_MPI,	size_y_linear_S_3_MPI, counts_G_lv_3_MPI, counts_G_lv_3_MPI_non_redundant ,counts_S_lv_3_MPI,		G3_linear_MPI,			y_linear_S_3_MPI,			lambda_log,			exposure_rate,			sigma_inv_log,			prop_3		),
					package_real(	ct_in_levels[4],	Q,	shards_in_levels[4],	size_counts_G_lv_4_MPI,	size_counts_G_lv_4_MPI_non_redundant,	size_counts_S_lv_4_MPI,	size_G4_linear_MPI,	size_y_linear_S_4_MPI, counts_G_lv_4_MPI, counts_G_lv_4_MPI_non_redundant, counts_S_lv_4_MPI,		G4_linear_MPI,			y_linear_S_4_MPI,			lambda_log,			exposure_rate,			sigma_inv_log,			prop_4		)
		),
		real_data2,
		int_package
	));
}
