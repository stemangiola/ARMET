functions{
	vector horseshoe_get_tp(vector zb, vector[] local, real[] global, real scale_global, real c2) {
		int K = rows(zb);
		vector[K] lambda = local[1] .* sqrt(local[2]);
		vector[K] lambda2 = square(lambda);
		real tau = global[1] * sqrt(global[2]) * scale_global;
		vector[K] lambda_tilde = sqrt(c2 * lambda2 ./ (c2 + tau^2 * lambda2));
		return zb .* lambda_tilde * tau;
	}

	real horseshoe_get_lp(vector zb, vector[] local, real df, real[] global, real df_global, real c2, real df_slab){

		// real<lower=0> hs_df; // == 1  // If divergencies increase this
		// real<lower=0> hs_df_global; // == 1
		// real<lower=0> hs_df_slab; // == 4 // df of the outliers

		vector[6] lp;

		lp[1] = normal_lpdf(zb | 0, 1);
		lp[2] = normal_lpdf(local[1] | 0, 1) - 101 * log(0.5);
		lp[3] = inv_gamma_lpdf(local[2] | 0.5 * df, 0.5 * df);
		lp[4] = normal_lpdf(global[1] | 0, 1)  - 1 * log(0.5);
		lp[5] = inv_gamma_lpdf(global[2] | 0.5 * df_global, 0.5 * df_global);
		lp[6] = inv_gamma_lpdf(c2 | 0.5 * df_slab, 0.5 * df_slab);

		return(sum(lp));
	}
}
data {
		// Horseshoe
	real<lower=0> hs_df; // == 1  // If divergencies increase this
	real<lower=0, upper=1> par_ratio; // real<lower=0> hs_scale_global; // Ratio of the expected number of non-zero coefficients     !! KEY PARAMETER
	real<lower=0> hs_scale_slab; // == 2 // regularisation/scale of outliers                !! KEY PARAMETER

}
	transformed data{
		int GM = 1000;
		real df_global = 3;
		real df_slab = 25;
	}
	parameters{

		vector<lower=0>[GM] error_ref_mix_z;

		// Horseshoe
		vector<lower=0>[GM] hs_local[2]; // local parameters for horseshoe prior
		real<lower=0> hs_global[2]; // horseshoe shrinkage parameters
		real<lower=0> hs_c2; // horseshoe shrinkage parameters
	}
	transformed parameters{
		// Horseshoe
		vector[GM] error_ref_mix = horseshoe_get_tp(error_ref_mix_z, hs_local, hs_global, par_ratio / sqrt(GM), hs_scale_slab^2 * hs_c2);


	}
model{
	target += horseshoe_get_lp(error_ref_mix_z, hs_local, hs_df, hs_global, df_global, hs_c2, df_slab);

}
