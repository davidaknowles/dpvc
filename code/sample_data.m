function Y=sample_data(param_set, N, D)

e= param_set.mu_noise+ param_set.sigma_noise*randn(D,N);
Y = (param_set.G.*param_set.C) * param_set.X + e;

end