function param_set = one_iteration(Y, param_set, settings)

param_set = sample_X(Y, param_set);
param_set = sample_C_collapsed(Y, param_set, settings);

if settings.sample_alpha
    param_set = slice_sampler_hyper_CRP(param_set);
end

if settings.sample_sigma_noise
    param_set = sample_sigma_noise(Y, param_set);
end

if settings.sample_sigma_g
    param_set = sample_sigma_g(Y,param_set);
end

end

