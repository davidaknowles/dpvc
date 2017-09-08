function [param] = sample_C_collapsed_d(Y, param, settings, d)
% this is just used by sequential initialisation
[D N] = size(Y);
K= size(param.C, 2);

Cs=[param.C; zeros(1,K)];
cc = param.cc;
Ns=param.Ns;
Gs=[param.G; param.mu_g+param.sigma_g * randn(1,K)];
Xs=param.X;

m_aux = settings.m_aux;

Xs = [Xs; param.mu_x+param.sigma_x * randn(m_aux, N)]; 
X2cache=sum(Xs.^2,2);

p= log(Ns.*(1/(param.alpha+D-1)));
p((end+1):(end+m_aux)) = log((param.alpha/m_aux)) -log(param.alpha+D-1);
Gsamps =zeros(1,K+m_aux);
xy=Xs*Y(d,:)'; 
% y2=Y(d,:)*Y(d,:)'; 
for k = 1:K+m_aux

    lambda = X2cache(k)/param.sigma_noise^2+1/param.sigma_g^2;
    mu=(xy(k)/param.sigma_noise^2+param.mu_g/param.sigma_g^2)/lambda; 
    p(k) = p(k) -.5*(log(lambda) - lambda*mu^2);
    % the following terms are constant wrt k
    % -.5*(N*log(2*pi)+N*log(param.sigma_noise^2) + Y2cache(d)/param.sigma_noise^2 ...
    %    + log(param.sigma_g^2) + param.mu_g^2/param.sigma_g^2 )
        
    Gsamps(k)=mu + lambda^(-.5) * randn(); 
end
p = exp(p-max(p));
p=p/(sum(p));

% sample from the conditional probabilities
uu = rand;
kk = 1+sum(uu>cumsum(p));
g=Gsamps(kk); 
if kk>K
    Xs(K+1,:)= Xs(kk,:);
    Cs=[Cs zeros(D,1)];
    kk = K+1;

    K=K+1;
    Ns(kk)=0;
end
Xs(K+1:end,:)=[]; 

Gs(d,kk)=g;
cc(d)=kk;
Ns(kk) = Ns(kk) + 1; % increment number of data items in component kk

Cs(d,:)=zeros(1,K);
Cs(d,kk)=1;
   
param.C=Cs;
param.cc=cc;
param.Ns=Ns;
param.G=Gs;
param.X=Xs;

end