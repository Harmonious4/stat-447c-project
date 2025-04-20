data {
  int N;
  vector[N] t;
  vector[N] y;
}

parameters {
  real intercept;
  real slope;
  real quadratic;
  real<lower=0,upper=1> previous_coef;
  real<lower=0> sigma;
  real sin_coef;
  real cos_coef;
}

transformed parameters {
  vector[N] delta = sin_coef*sin(2*pi()*t) + cos_coef*cos(2*pi()*t);
  vector[N] mu = intercept + slope*t + quadratic*t^2 + delta;
}

model {
  intercept ~ normal(10, 10);
  slope ~ normal(0, 0.1);
  quadratic ~ normal(0, 0.01);
  previous_coef ~ beta(2, 2);
  sigma ~ exponential(0.5);
  sin_coef ~ normal(0, 10);
  cos_coef ~ normal(0, 10);
  y[1] ~ normal(mu[1], sigma);
  for (i in 2:N) {
    y[i] ~ normal(mu[i] + previous_coef^((t[i]-t[i-1])*365.25)*(y[i-1]-mu[i-1]), sigma);
  }
}

generated quantities {
  real estimate2025 = intercept + slope * 25 + quadratic * 25^2;
  real slope2025 = slope + 25 * quadratic;
  real amplitude = hypot(sin_coef, cos_coef);
  real min_day = atan2(-sin_coef, amplitude - cos_coef) * 365.25 / pi();
}
