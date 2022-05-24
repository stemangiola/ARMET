functions{
	vector[] dirichlet_regression_rng( matrix X, matrix beta, real phi, real plateau){
		
		vector[cols(beta)] p[rows(X)];
		matrix[num_elements(p[1]), num_elements(p[,1])] mu = (X * beta)';
    real buffer;
    
		for(i in 1:cols(mu)) {
			mu[,i] = softmax(mu[,i]);
			buffer = 1.0/min(mu[,i]) * plateau;
			
			p[i] = dirichlet_rng(mu[,i] * phi * buffer);
			
		}
		
		return (p);
		
	}
	

  matrix get_mean_prop(matrix X, matrix alpha){

	  	matrix[rows(X), cols(alpha)] mu = X * alpha;

		  for(j in 1:num_elements(mu[,1])) mu[j]  = to_row_vector(softmax( to_vector( mu[j] )));
  
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


	// lv4
  matrix[A * (lv == 4) ,ct_in_nodes[8]]  alpha_g; // dendritic myeloid
  matrix[A * (lv == 4) ,ct_in_nodes[9]]  alpha_h; // macrophage
  matrix[A * (lv == 4) ,ct_in_nodes[10]]  alpha_i; // NK
  matrix[A * (lv == 4) ,ct_in_nodes[11]] alpha_l; // CD4
  matrix[A * (lv == 4) ,ct_in_nodes[12]] alpha_m; // CD8

	real<lower=0> phi[12]; //[fam_dirichlet ? 10 : ct_in_levels[lv]];

	matrix[Q,A] X_scaled;
}
generated quantities{

	// lv4
  vector[ct_in_nodes[8]]  prop_g_rng[Q * (lv == 4)]; // dendritic myeloid childrens
  vector[ct_in_nodes[9]]  prop_h_rng[Q * (lv == 4)]; // macrophage childrens
  vector[ct_in_nodes[10]] prop_i_rng[Q * (lv == 4)]; // nk primed
  vector[ct_in_nodes[11]] prop_l_rng[Q * (lv == 4)]; // CD4 childrens
  vector[ct_in_nodes[12]] prop_m_rng[Q * (lv == 4)]; // CD8 childrens
  
  matrix[Q * (lv == 4), ct_in_nodes[8]]  mu_g_rng; // dendritic myeloid childrens
  matrix[Q * (lv == 4),ct_in_nodes[9]]  mu_h_rng; // macrophage childrens
  matrix[Q * (lv == 4),ct_in_nodes[10]] mu_i_rng; // nk primed
  matrix[Q * (lv == 4),ct_in_nodes[11]] mu_l_rng; // CD4 childrens
  matrix[Q * (lv == 4),ct_in_nodes[12]] mu_m_rng; // CD8 childrens
  

	prop_g_rng = dirichlet_regression_rng( X_scaled, alpha_g, (phi[1]) , 1);
	prop_h_rng = dirichlet_regression_rng( X_scaled, alpha_h, (phi[2]) , 1);
	prop_i_rng = dirichlet_regression_rng( X_scaled, alpha_i, (phi[3]) , 1);
	prop_l_rng = dirichlet_regression_rng( X_scaled, alpha_l, (phi[4]) , 1);
	prop_m_rng = dirichlet_regression_rng( X_scaled, alpha_m, (phi[5]) , 1);

	mu_g_rng = get_mean_prop(X_scaled, alpha_g);
	mu_h_rng = get_mean_prop(X_scaled, alpha_h);
	mu_i_rng = get_mean_prop(X_scaled, alpha_i);
	mu_l_rng = get_mean_prop(X_scaled, alpha_l);
	mu_m_rng = get_mean_prop(X_scaled, alpha_m);

}
