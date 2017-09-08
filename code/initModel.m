function param=initModel(settings, param)
% Make an initial sample for the MCMC sampler.
if settings.sample_alpha
    param.alpha=gamrnd(1,1);
end

if settings.sample_sigma_g
    variance=1/gamrnd(param.sigma_g_hyper_a, 1/param.sigma_g_hyper_b );
    param.sigma_g=sqrt(variance);
end

if settings.K
    [param.C param.cc param.Ns]= crpgen(settings.D, param.alpha, settings.K);
else
    [param.C param.cc param.Ns]= crpgen(settings.D, param.alpha);
    K=size(param.C,2);
    assert(sum(sum(param.C))==settings.D);
end

param.G = param.mu_g+param.sigma_g*randn(settings.D , K);
param.X = param.mu_x + param.sigma_x.*randn(K, settings.N);

param.mu_noise=0;

if settings.sample_sigma_noise
    variance=1/gamrnd(param.noise_hyper_a,1/param.noise_hyper_b );
    param.sigma_noise=sqrt(variance);
end

end

