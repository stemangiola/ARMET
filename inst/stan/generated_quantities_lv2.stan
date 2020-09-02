functions{
			vector dirichlet_regression_rng( row_vector X, matrix alpha, real phi, real plateau){

		// Calculate log prob
		return (dirichlet_rng( softmax( append_row([0]', to_vector(X * alpha))) * exp( phi) + plateau ));
}


vector[] beta_regression_rng( matrix X, matrix alpha, real[] phi, real plateau){

		vector[cols(alpha)+1] p[rows(X)];

		//matrix[num_elements(p[,1]), num_elements(p[1])] mu;
		vector[num_elements(phi)] phi_exp= ( 1.0 ./ (to_vector(phi )+ 0.0001));

		for(j in 1:num_elements(p[,1])) {

			vector[num_elements(p[1])] mu  = softmax( append_row([0]', to_vector(X[j] * alpha)));
		
		p[j] = to_vector(beta_rng((mu .* phi_exp) +plateau, ((1.0 - mu) .* phi_exp) + plateau));
		
	}
	return (p);
}

  vector[] get_mean_prop(matrix X, matrix alpha){

	  	vector[cols(alpha)+1] mu[rows(X)];

		  for(j in 1:num_elements(mu[,1])) 
	  		mu[j]  = softmax( append_row([0]', to_vector(X[j] * alpha)));
  
  	return(mu);
  }

} 
data {
	// shards

	int lv;
	// Reference matrix inference

	// Cell types
 	int<lower=0> Q;
  int<lower=1> n_nodes;
  int<lower=1> ct_in_nodes[n_nodes];

	// Dirichlet regression
	int A; // factors of interest

} parameters {



  matrix[A,ct_in_nodes[2]-1]  alpha_a; // Immune cells


	real<lower=0> phi[12]; //[fam_dirichlet ? 10 : ct_in_levels[lv]];
	matrix[Q,A] X_scaled;




}
generated quantities{


  vector[ct_in_nodes[2]]  prop_a_rng[Q * (lv == 2)]; // Immune cells childrens
  vector[ct_in_nodes[2]]  mu_a_rng[Q * (lv == 2) ]; 


  for(q in 1:Q) prop_a_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_a, phi[1] , 0.05);
  mu_a_rng = get_mean_prop(X_scaled, alpha_a);


}
