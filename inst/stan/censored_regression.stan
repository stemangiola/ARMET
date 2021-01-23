functions {
		real censored_regression_single_lpdf(real time, vector prop, matrix alpha, vector phi, int is_censored){
		// NOT USED
		int C = rows(prop);
		vector[C] mu;
		real lp = 0;
	
		// For each cell type
		for(c in 1:C) mu[c] = alpha[1,c] + (prop[c]) * alpha[2,c];
	
	  // special treatment of censored data
	  if (is_censored == 0) lp += normal_lpdf(rep_vector(time, C) | mu, phi);
	  else if (is_censored == 1) lp += normal_lccdf(rep_vector(time, C) | mu, phi); 
		else print("Wrong cens value!!!");
		
		return(lp);
	  
	}
	
	real censored_regression_lpdf(vector time, matrix[] X, matrix alpha, row_vector phi, int[] which_censored, int[] which_non_censored){
		
		int C = size(X);
		int S = rows(X[1]);
		matrix[S,C] mu;
		real lp = 0;
	
		// For each cell type
	for(c in 1:C) mu[,c] = X[c] * alpha[,c];
	
	 // print((rep_matrix(time[which_censored], C)));
	 // print((mu[which_censored]));
	 // print((rep_matrix(phi, size(which_censored))));
	 // print("---");
	
	// special treatment of censored data
  lp += normal_lpdf(to_vector(rep_matrix(time[which_non_censored], C)) | to_vector((mu[which_non_censored])), to_vector((rep_matrix(phi, size(which_non_censored)))));
  lp += normal_lccdf(to_vector(rep_matrix(time[which_censored], C)) | to_vector((mu[which_censored])), to_vector((rep_matrix(phi, size(which_censored)))));

	return(lp);
	  
	}
	
}
data {
  // Cell types
 	int<lower=0> S;
 	int<lower=1> C;
 	int A; // factors of interest
 	
	vector[S] time;
	matrix[S, A] X[C];
	
	
	// Censoring
	int	n_cens;
	int n_non_cens;
	int	which_censored[n_cens];
	int which_non_censored[n_non_cens];

}
parameters {
  matrix[A ,C]  alpha;
  row_vector<lower=0>[C] phi;
}
model {
	// print(phi);
	// print(alpha);
	
  time ~ censored_regression(X, alpha, phi, which_censored, which_non_censored);
  to_vector(alpha) ~ student_t(3, 0, 5);
  phi ~ normal(0,3); //student_t(3, 0, 5);
  
}

