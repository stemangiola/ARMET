data {
  int<lower=0> N;
  int<lower=0> counts[N];
}

parameters {
  real<lower=0> mu;
  real<lower=0> phi;
}

model {
  counts ~ neg_binomial_2(mu, phi);
}

