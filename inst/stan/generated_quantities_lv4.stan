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


	// lv4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[8]-1]  alpha_g; // dendritic myeloid
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[9]-1]  alpha_h; // macrophage
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[10]-1]  alpha_i; // NK
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[11]-1] alpha_l; // CD4
  matrix[A * (lv == 4) * do_regression,ct_in_nodes[12]-1] alpha_m; // CD8

	real<lower=0> phi[12]; //[fam_dirichlet ? 10 : ct_in_levels[lv]];

	matrix[Q,A] X_scaled;
}
generated quantities{

	// lv4
  vector[ct_in_nodes[8]]  prop_g_rng[Q * (lv == 4)* do_regression]; // dendritic myeloid childrens
  vector[ct_in_nodes[9]]  prop_h_rng[Q * (lv == 4)* do_regression]; // macrophage childrens
  vector[ct_in_nodes[10]] prop_i_rng[Q * (lv == 4)* do_regression]; // nk primed
  vector[ct_in_nodes[11]] prop_l_rng[Q * (lv == 4)* do_regression]; // CD4 childrens
  vector[ct_in_nodes[12]] prop_m_rng[Q * (lv == 4)* do_regression]; // CD8 childrens
  
  vector[ct_in_nodes[8]]  mu_g_rng[Q * (lv == 4)* do_regression]; // dendritic myeloid childrens
  vector[ct_in_nodes[9]]  mu_h_rng[Q * (lv == 4)* do_regression]; // macrophage childrens
  vector[ct_in_nodes[10]] mu_i_rng[Q * (lv == 4)* do_regression]; // nk primed
  vector[ct_in_nodes[11]] mu_l_rng[Q * (lv == 4)* do_regression]; // CD4 childrens
  vector[ct_in_nodes[12]] mu_m_rng[Q * (lv == 4)* do_regression]; // CD8 childrens
  

	for(q in 1:Q) prop_g_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_g, phi[1] , 1);
	for(q in 1:Q) prop_h_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_h, phi[2] , 1);
	for(q in 1:Q) prop_i_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_i, phi[3] , 1);
	for(q in 1:Q) prop_l_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_l, phi[4] , 1);
	for(q in 1:Q) prop_m_rng[q] = dirichlet_regression_rng( X_scaled[q], alpha_m, phi[5] , 1);

	mu_g_rng = get_mean_prop(X_scaled, alpha_g);
	mu_h_rng = get_mean_prop(X_scaled, alpha_h);
	mu_i_rng = get_mean_prop(X_scaled, alpha_i);
	mu_l_rng = get_mean_prop(X_scaled, alpha_l);
	mu_m_rng = get_mean_prop(X_scaled, alpha_m);

}
