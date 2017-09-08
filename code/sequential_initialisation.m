function param = sequential_initialisation(Y, param, settings)

[D N]=size(Y); 
param.X=randn(1,N);
param.G=randn(1,1); 
param.C=ones(1,1); 
param.cc=[1]; 
param.Ns=[1]; 

Yd=Y(1,:); 
true_mask = settings.mask_unObs;
settings.mask_unObs = true_mask(1,:); 

Yd(settings.mask_unObs)=0; 
param = one_iteration(Yd, param, settings);

for d = 2:D
    Yd=[Yd; Y(d,:)];
    Yd(end,true_mask(d,:))=0.0; 
    settings.mask_unObs = true_mask(1:d,:); 

    % sample the unobserved data from the model
    [param] = sample_C_collapsed_d(Yd, param, settings, d);
    if settings.sequential_expensive
    	param = one_iteration(Yd, param, settings);
    end
    Yd= sample_unObs(Yd, param, settings);

    size(param.C,2);
end