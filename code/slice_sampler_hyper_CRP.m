function [ param ] = slice_sampler_hyper_CRP(param)

alpha=param.alpha;
Cs=param.C;
Ns=param.Ns;

%p(alpha | z) prop_to p(z|alpha)*p(alpha)
%p(alpha)=exp(-alpha) --> log_prior=-alpha

[D K]=size(Cs);
temp = sum(gammaln(Ns));

f= @(g) gammaln(g)+K*log(g) +temp - gammaln(g+D);
h =@(x) f(exp(x))+x-exp(x);
N = 1;

t= slicesample(log(alpha),N,'logpdf',h,'thin',1,'burnin',1);
param.alpha=exp(t);

end