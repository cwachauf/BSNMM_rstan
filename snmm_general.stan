data
{
  int<lower=1> N; // number of data points
  int<lower=1> M; // number of mixtures/components
  real y[N]; // data points
  vector<lower=0>[M] alphas; // priors for the mixture components
  real<lower=0> mus0[M]; // priors for the mixture means
  real<lower=0> std0; // the hyperparameter: standard-deviation for the means...
}

parameters
{
  simplex[M] phis; 
  real locations[M];
  real<lower=0> scales[M];
  real shapes[M];
}

model
{
  real ps[M];
  locations ~ normal(mus0,std0);
  scales ~ cauchy(0,2.5);
  shapes ~ normal(0,0.5);
  phis ~ dirichlet(alphas);
  for(n in 1:N)
  {
    for(i in 1:M)
    {
      ps[i] <- log(phis[i]) + skew_normal_log(y[n],locations[i],scales[i],shapes[i]);
    }
    increment_log_prob(log_sum_exp(ps));
  }
}