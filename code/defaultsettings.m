function [settings param]=defaultsettings()

settings.sequential_expensive=1; 
settings.save_all=0; 
settings.iterations=100;
settings.burnin=0; 
settings.thinout=1; 
settings.m_aux=3;
settings.K=0;

settings.sample_alpha = 0;
settings.sample_sigma_x=0;
settings.sample_sigma_g=0;
settings.sample_sigma_noise=0;

settings.output_file='output.mat';

param.alpha=1;
param.mu_x=0;
param.sigma_x=1;
param.mu_g=0;
param.sigma_g=1;
param.mu_noise=0;
param.sigma_noise=0.1;
param.noise_hyper_a=1;
param.noise_hyper_b=.1;
param.sigma_g_hyper_a=1;
param.sigma_g_hyper_b=1;

end
