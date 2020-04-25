functions{
  
  real beta_regression_logit_lpdf(vector[] p, matrix X, matrix alpha, matrix beta){

		real lp = 0;
		matrix[num_elements(p[,1]), num_elements(p[1])] mu;
		matrix[num_elements(p[,1]), num_elements(p[1])] phi;

		mu  = inv_logit(X * alpha);
		phi  = exp(X * beta);

		for(j in 1:num_elements(p[1])) {

      	 lp += beta_lpdf(p[,j] | (mu[,j] .* phi[,j]) , ((1.0 - mu[,j]) .* phi[,j]) );

		}
		return (lp);
	}
	
	real beta_regression_softmax_lpdf(vector[] p, matrix X, matrix alpha, matrix beta){

		real lp = 0;
		
		for(j in 1:num_elements(p[,1])) {

			vector[num_elements(p[1])] mu  = softmax( append_row([0]', to_vector(X[j] * alpha)));
			vector[num_elements(p[1])] phi  = exp( to_vector(X[j] * beta));

     	lp += beta_lpdf(p[j] | (mu .* phi), ((1.0 - mu) .* phi) );

		}
		return (lp);
	}
	
}
// The input data is a vector 'y' of length 'N'.
data {
  
  int S;
  int A; // factors of interest
  int C;
	matrix[S,A] X;
  simplex[C]  prop[S];
  
}

// The parameters accepted by the model. Our model
// accepts two parameters 'mu' and 'sigma'.
parameters {
  matrix[A,C-1] alpha_generative;
  matrix[A,C] alpha_descriptive;
  matrix[A,C] beta_descriptive;
  matrix[A,C] beta_generative;
  
}


// The model to be estimated. We model the output
// 'y' to be normally distributed with mean 'mu'
// and standard deviation 'sigma'.
model {
  to_vector(beta_descriptive) ~ normal(0, 10);
  to_vector(beta_generative) ~ normal(0, 10);
  to_vector(alpha_descriptive) ~ normal(0, 10);
  to_vector(alpha_generative) ~ normal(0, 10);

  prop ~ beta_regression_logit(X, alpha_descriptive, beta_descriptive);
  prop ~ beta_regression_softmax(X, alpha_generative, beta_generative);

}

