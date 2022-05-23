functions{

vector[] beta_regression_rng( matrix X, matrix beta, vector phi){

 vector[cols(beta)] p[rows(X)];
 matrix[num_elements(p[1]), num_elements(p[,1])] mu = (X * beta)';

 
    real buffer;
    
		for(i in 1:cols(mu)) mu[,i] = softmax(mu[,i]);
		for(i in 1:cols(mu)) {
			for(j in 1:rows(mu)) {
			
			buffer = 1.0/fmin(mu[j,i], 1-mu[j,i]) * 0.5;
			
			p[i,j] = beta_rng(mu[j,i] .* phi[j] * buffer, (1.0 - mu[j,i]) .* phi[j] * buffer);
			
		}}
		
		return (p);
	}


  vector[] get_mean_prop(matrix X, matrix alpha){

	  	vector[cols(alpha)] mu[rows(X)];

		  for(j in 1:num_elements(mu[,1])) 
	  		mu[j]  = softmax( to_vector(X[j] * alpha));
  
  	return(mu);
  }

} 
data {
	// shards

	// Reference matrix inference

	// Cell types
 	int<lower=0> Q;
  int<lower=2> number_of_cell_types;

	// Dirichlet regression
	int A; // factors of interest

} 
parameters {

  matrix[A,number_of_cell_types]  alpha_1; // Root
	vector<lower=0>[number_of_cell_types] phi; //[fam_dirichlet ? 10 : ct_in_levels[lv]];
	matrix[Q,A] X_scaled;

}
generated quantities{

  // lv1
  vector[number_of_cell_types]  prop_1_rng[Q  ]; // Root
  vector[number_of_cell_types]  mu_1_rng[Q  ]; 



	prop_1_rng = beta_regression_rng(X_scaled, alpha_1, phi	);
	mu_1_rng = get_mean_prop(X_scaled, alpha_1);

 
}
