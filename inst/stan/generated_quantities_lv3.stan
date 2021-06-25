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



  // lv3
  matrix[A,ct_in_nodes[3]]  alpha_b; // b cells
  matrix[A,ct_in_nodes[4]]  alpha_c; // granulocyte
  matrix[A,ct_in_nodes[5]]  alpha_d; // mono_derived
  matrix[A,ct_in_nodes[6]]  alpha_e; // natural_killer
  matrix[A,ct_in_nodes[7]]  alpha_f; // t_cell


	vector<lower=0>[12] phi; //[fam_dirichlet ? 10 : ct_in_levels[lv]];
	matrix[Q,A] X_scaled;




}
generated quantities{


  // lv3
  vector[ct_in_nodes[3]]  prop_b_rng[Q * (lv == 3)]; // b cells childrens
  vector[ct_in_nodes[4]]  prop_c_rng[Q * (lv == 3)]; // granulocyte childrens
  vector[ct_in_nodes[5]]  prop_d_rng[Q * (lv == 3)]; // mono_derived childrens
  vector[ct_in_nodes[6]]  prop_e_rng[Q * (lv == 3)]; // natural_killer childrens
  vector[ct_in_nodes[7]]  prop_f_rng[Q * (lv == 3)]; // t_cell childrens

  matrix[Q * (lv == 3) ,ct_in_nodes[3]] mu_b_rng; 
  matrix[Q * (lv == 3) ,ct_in_nodes[4]]  mu_c_rng; 
  matrix[Q * (lv == 3) ,ct_in_nodes[5]]  mu_d_rng; 
  matrix[Q * (lv == 3) ,ct_in_nodes[6]]  mu_e_rng; 
  matrix[Q * (lv == 3) ,ct_in_nodes[7]]  mu_f_rng; 


	prop_b_rng = dirichlet_regression_rng( X_scaled, alpha_b, exp(phi[1]) , 1);
	prop_c_rng = dirichlet_regression_rng( X_scaled, alpha_c, exp(phi[2]) , 1);
	prop_d_rng = dirichlet_regression_rng( X_scaled, alpha_d, exp(phi[3]) , 1);
	prop_e_rng = dirichlet_regression_rng( X_scaled, alpha_e, exp(phi[4]) , 1);
	prop_f_rng = dirichlet_regression_rng( X_scaled, alpha_f, exp(phi[5]) , 1);

	mu_b_rng = get_mean_prop(X_scaled, alpha_b);
	mu_c_rng = get_mean_prop(X_scaled, alpha_c);
	mu_d_rng = get_mean_prop(X_scaled, alpha_d);
	mu_e_rng = get_mean_prop(X_scaled, alpha_e);
	mu_f_rng = get_mean_prop(X_scaled, alpha_f);

}
