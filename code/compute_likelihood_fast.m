function [ lh ] = compute_likelihood_fast(Y,  param)

[D N] = size(Y);
Gm=param.G.*param.C;
E=Y-Gm*param.X;
lh=-D*N*(log(param.sigma_noise)+0.5*log(2.0*pi)) - 0.5*sum(sum(E.^2))/(param.sigma_noise^2); 

end