function [param] = sample_sigma_g(Y, param)
G=param.G;
C=param.C;
Gm=G.*C;

[D N] = size(Y); 
variance=1/gamrnd(param.sigma_g_hyper_a+D/2, 1/(param.sigma_g_hyper_b+.5*sum(sum(Gm.*Gm)))); 
param.sigma_g=sqrt(variance);

end