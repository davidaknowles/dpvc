function [energy] = compute_energy(Y, param)

K= size(param.C,2);
[D N] = size(Y);

lk =  compute_likelihood_fast(Y, param);

prior_g= -0.5*(1/(param.sigma_g^2)) *(sum(sum( (param.G-repmat(param.mu_g.*ones(D,1),[1 K])).^2  ) ));
prior_x = -0.5*(1/(param.sigma_x^2)) *(sum(sum( (param.X-repmat(param.mu_x.*ones(K,1),[1 N])).^2  ) ));
prior_c = gammaln(param.alpha)+K*log(param.alpha) +sum(gammaln(param.Ns)) - gammaln(param.alpha+D);

energy = lk+prior_g + prior_x +prior_c;

end