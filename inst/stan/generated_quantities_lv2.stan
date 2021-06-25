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



  matrix[A,ct_in_nodes[2]]  alpha_a; // Immune cells


	vector<lower=0>[12] phi; //[fam_dirichlet ? 10 : ct_in_levels[lv]];
	matrix[Q,A] X_scaled;




}
generated quantities{


  vector[ct_in_nodes[2]]  prop_a_rng[Q * (lv == 2)]; // Immune cells childrens
  matrix[Q * (lv == 2), ct_in_nodes[2]]  mu_a_rng; 


  prop_a_rng = dirichlet_regression_rng( X_scaled, alpha_a, exp(phi[1]) , 0.5);
  mu_a_rng = get_mean_prop(X_scaled, alpha_a);


}
