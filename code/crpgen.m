function [Cs cc Ns] = crpgen(D, alpha, maxK)

if nargin==2
    maxK=Inf;
end

K=0;
Ns=zeros(1,K);
cc=zeros(1,D);
for d=1:D
    p= log(Ns.*(1/(alpha+D-1)));

    if K<maxK
        p(end+1) = log(alpha) -log(alpha+D-1);
    end

    p = exp(p-max(p));
    p=p/(sum(p));

    uu = rand;
    kk = 1+sum(uu>cumsum(p));

    if kk>K
        kk = K+1;
        K=K+1;
        Ns(kk)=0;
    end

    cc(d)=kk;
    Ns(kk) = Ns(kk) + 1;

end

Cs=zeros(D,K);
cc_ind = sub2ind(size(Cs), 1:D, cc );
Cs(cc_ind)=1;

end