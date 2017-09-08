function [param] = sample_sigma_noise(Y, param)
G=param.G;
C=param.C;
X=param.X;
[D N] = size(Y); 
E = Y - (G.*C)*X;
% trace(E'*E) = sum(sum(E.*E));
variance=1./gamrnd((param.noise_hyper_a)+(D*N)/2,1/(param.noise_hyper_b+.5*trace(E'*E)));
param.sigma_noise=sqrt(variance);

end