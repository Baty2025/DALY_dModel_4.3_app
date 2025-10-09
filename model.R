library(rstan)

stan_model_tp <- stan_model(model_code = "
data {
  int<lower=0> T_pos;
  int<lower=0> N;
  real<lower=0, upper=1> Se;
  real<lower=0, upper=1> Sp;
}
parameters {
  real<lower=0, upper=1> pi;
}
model {
  T_pos ~ binomial(N, pi * Se + (1 - pi) * (1 - Sp));
  pi ~ beta(1, 1);
}
")

saveRDS(stan_model_tp, file = "stan_model_tp.rds")
