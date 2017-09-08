% This script demoes the MCMC sampler for the Dirichlet Process Variable 
% Clustering model described in 
%
% Konstantina Palla, David A. Knowles and Zoubin Ghahramani
% A nonparametric variable clustering model
% NIPS 2012
% 
% NOTE: if you run this on your own data we would recommend pre-normalising
% i.e. center and scale
rng(12345); 

% generate synthetic data (the first 20 dimensions are in one cluster,
% the next 20 are in another)
D=40;
K=2;
N=100;
C = zeros(D,K);
C(1:20,1)=1;
C(21:40,2)=1;
cc=zeros(1,D);
cc(1:20)=1;
cc(21:40)=2;

[settings param]=defaultsettings();
settings.sample_alpha = 1;
settings.sample_sigma_x=1;
settings.sample_sigma_g=1;
settings.sample_sigma_noise=1;
settings.iterations=1000;

G=param.mu_g+param.sigma_g*randn(D,K);
X=param.mu_x + param.sigma_x*randn(K,N);

e=param.mu_noise+param.sigma_noise*randn(D,N);

Y = (G.*C)*X+e;

[settings.D settings.N]=size(Y);

% use some heldout data
hold_out_per=0.1;
holdout=ceil(hold_out_per*(settings.D*settings.N));
per=randperm(settings.D*settings.N);
mask_unObs = logical(sparse( zeros(size(Y))));
mask_unObs(per(1:holdout))=1;
settings.mask_unObs=mask_unObs;

% run inferenec
tic
[ param ] = sequential_initialisation(Y, param,  settings);
toc

tic
[ param ] = inference(Y, param,  settings);
toc

% plot the energy over time
plot(param.E_iter)

% print out the best clustering found
disp(param.bestSample.cc)