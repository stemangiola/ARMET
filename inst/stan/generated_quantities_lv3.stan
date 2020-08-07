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

	int<lower=0> S;
	int CL; // counts linear size

	// Cell types
 	int<lower=0> Q;
  int<lower=1> n_nodes;
  int<lower=1> ct_in_nodes[n_nodes];

	// Dirichlet regression
	int A; // factors of interest
	matrix[Q,A] X;
	int do_regression;

} parameters {



  // lv3
  matrix[A,ct_in_nodes[3]-1]  alpha_b; // b cells
  matrix[A,ct_in_nodes[4]-1]  alpha_c; // granulocyte
  matrix[A,ct_in_nodes[5]-1]  alpha_d; // mono_derived
  matrix[A,ct_in_nodes[6]-1]  alpha_e; // natural_killer
  matrix[A,ct_in_nodes[7]-1]  alpha_f; // t_cell


	real<lower=0> phi[12]; //[fam_dirichlet ? 10 : ct_in_levels[lv]];
	matrix[Q,A] X_scaled;




}
generated quantities{


  // lv3
  vector[ct_in_nodes[3]]  prop_b_rng[Q * (lv == 3)* do_regression]; // b cells childrens
  vector[ct_in_nodes[4]]  prop_c_rng[Q * (lv == 3)* do_regression]; // granulocyte childrens
  vector[ct_in_nodes[5]]  prop_d_rng[Q * (lv == 3)* do_regression]; // mono_derived childrens
  vector[ct_in_nodes[6]]  prop_e_rng[Q * (lv == 3)* do_regression]; // natural_killer childrens
  vector[ct_in_nodes[7]]  prop_f_rng[Q * (lv == 3)* do_regression]; // t_cell childrens

  vector[ct_in_nodes[3]]  mu_b_rng[Q * (lv == 3) * do_regression]; 
  vector[ct_in_nodes[4]]  mu_c_rng[Q * (lv == 3) * do_regression]; 
  vector[ct_in_nodes[5]]  mu_d_rng[Q * (lv == 3) * do_regression]; 
  vector[ct_in_nodes[6]]  mu_e_rng[Q * (lv == 3) * do_regression]; 
  vector[ct_in_nodes[7]]  mu_f_rng[Q * (lv == 3) * do_regression]; 


	for(q in 1:Q) prop_b_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_b, phi[1] , 1);
	for(q in 1:Q) prop_c_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_c, phi[2] , 1);
	for(q in 1:Q) prop_d_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_d, phi[3] , 1);
	for(q in 1:Q) prop_e_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_e, phi[4] , 1);
	for(q in 1:Q) prop_f_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_f, phi[5] , 1);

	mu_b_rng = get_mean_prop(X_scaled, alpha_b);
	mu_c_rng = get_mean_prop(X_scaled, alpha_c);
	mu_d_rng = get_mean_prop(X_scaled, alpha_d);
	mu_e_rng = get_mean_prop(X_scaled, alpha_e);
	mu_f_rng = get_mean_prop(X_scaled, alpha_f);

}
