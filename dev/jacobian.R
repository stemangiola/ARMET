library(rstan)

m_rate_wrong = "
transformed data{
	vector[4] alpha = [1,2,3,4]';
}
parameters{
	vector[3] rate_raw;
	
}
transformed parameters{
	vector[4] rate;

	simplex[4] prop;
	
	rate[1:3] = rate_raw;
	rate[4] = -sum(rate[1:3]);
	
	prop = softmax(rate);
}
model{
	prop ~ dirichlet(alpha);
}
"

fit_rate_wrong = stan(model_code = m_rate_wrong)

# Prop right
m_prop_right = "
transformed data{    vector[4] alpha = [1,2,3,4]';     }
parameters{	         simplex[4] prop;                  }
model{               prop ~ dirichlet(alpha);          }
"

fit_prop_right = stan(model_code = m_prop_right, cores=4)




# Jacobian
m_rate_right = "
functions{
vector sum_to_zero(vector v){
	int K = rows(v)+1;
	vector[K] v_0;
	
	v_0[1:(K-1)] = v;
	v_0[K] = -sum(v_0[1:(K-1)]);
	
	return(v_0);
}
}
transformed data{
	int K = 4;
	int N = 1;

	vector[K] alpha = [1,2,3,4]';
	matrix[K,K] ident = diag_matrix(rep_vector(1,K));

}
parameters{
	vector[K-1] rate_raw;
}
transformed parameters{
	simplex[K] prop[N];
	prop[1] = softmax(sum_to_zero(rate_raw));
}
model{

for(n in 1:N) {
	matrix[K,K] jacobian;

	prop[n] ~ dirichlet(alpha);
	
	for(i in 1:K) {
		for(j in 1:K) {
			jacobian[i,j] = prop[n,i] * (ident[i,j] - prop[n,j]);
		}
	}
	target += log_determinant(jacobian);
}

}
"

fit_rate_right = stan(model_code = m_rate_right, cores=4)
