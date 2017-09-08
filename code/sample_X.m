function [param] = sample_X(Y, param)

Gs=param.G;
Cs=param.C;

[D N] = size(Y);
K= size(Cs, 2);

Gm = Gs.*Cs;
tem=Gm' / param.sigma_noise^2;
lambda_diag= ((param.sigma_noise)^(-2))*sum(Gm.*Gm,1)+param.sigma_x^(-2);

for n=1:N
    mu = (tem*Y(:,n))./lambda_diag';%Dx1
    param.X(:,n) =  mu + lambda_diag'.^(-0.5).*randn(K,1);
end

% non-optimised version O(K^3) ------------------------------------------
% lambdae=ones(D,1)*param.sigma_noise^-2;
% lambdaG=bsxfun(@times,Gm,lambdae);
% prec=Gm'*lambdaG+eye(K) * param.sigma_x^2;
% mus=prec \ lambdaG'*Y;
% param.X=mus+chol(prec) \ randn(K,N);
