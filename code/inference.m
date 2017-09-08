
function [ param_cell] = inference(Y, param, settings)

iterations = settings.iterations;

%%%%%%%%%%%%%%%%%%%%%%%%%%% Initialize Markov Chain State %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% structure to store the Cs matrices over the iterations
cc_iter = cell(iterations , 1); %vector 1xD
Cs_iter = cell(iterations , 1); %matrix
assert(length(unique(cc_iter{1})) == size(Cs_iter{1},2));

%structure to store the number of dimensions assigned to each cluster
% over the iterations
Ns_iter  = cell(iterations, 1);
Gs_iter=cell(iterations,1);
Xs_iter=cell(iterations,1);
a_iter= zeros(iterations,1);
sigma_noise_iter=zeros(iterations,1);
sigma_g_iter=zeros(iterations,1);

% structure to store the energy
E_iter = zeros(iterations,1);

logpred = zeros(iterations,1);

% STORE THE real values of the Y_unObs
Y_true = Y;

tic

for i=1:settings.burnin
    Y= sample_unObs(Y, param, settings);
    param = one_iteration(Y, param, settings);
    fprintf('burnin it: %d K: %d noise std: %d g std: %d\n', i, size(param.C,2), param.sigma_noise, param.sigma_g);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Run the sampler %%%%%%%%%%%%%%%%%%%%%%%% 
[E_iter(1)] = compute_energy(Y, param);

Ns_iter{1}=param.Ns;
Gs_iter{1} = param.G;
Xs_iter{1} = param.X;
Cs_iter{1}=param.C;
cc_iter{1}=param.cc;
a_iter(1) = param.alpha;
sigma_noise_iter(1) = param.sigma_noise;
sigma_g_iter(1) = param.sigma_g;

test_mask = settings.mask_unObs & ~isnan(Y_true); 
Ntest = sum(sum(test_mask)); 

bestSample=param; 
bestE = E_iter(1); 

for i = 2:iterations

    for j=1:settings.thinout
        % sample the unobserved data from the model
        Y= sample_unObs(Y, param, settings);
        param = one_iteration(Y, param, settings);
    end
    
    Xs_iter{i}=param.X;
    Gs_iter{i}=param.G;
    Cs_iter{i}=param.C;
    cc_iter{i}=param.cc;
    Ns_iter{i}=param.Ns;
    a_iter(i)=param.alpha;
    sigma_noise_iter(i)=param.sigma_noise;
    sigma_g_iter(i)=param.sigma_g;
    timing(i)=toc(); 
    
    %compute energy function
    [E_iter(i)] = compute_energy(Y, param);
    
    if E_iter(i) > bestE
        bestE = E_iter(i); 
        bestSample = param; 
    end
     
    if (Ntest > 0)
        Ypred=(param.G.*param.C) * param.X;
        er = Ypred - Y_true; 
        logpred(i) = - .5 * ( Ntest * log(2.0 * pi * param.sigma_noise^2) + sum(sum( er(test_mask) .^ 2) )/ param.sigma_noise^2 ) ;
        logpred(i) = logpred(i); 
    end
    
%     plot_matrices(Cs_iter{i}, i);
    if (mod(i,1)==0)

        fprintf('it: %d (thin %d) time %d energy: %d lp: %d K: %d noise std: %d g std: %d\n', i, settings.thinout, timing(i), E_iter(i), logpred(i) / Ntest, size(param.C,2), param.sigma_noise, param.sigma_g);
        if settings.save_all
            save(settings.output_file); 
        else
            save(settings.output_file,'logpred','timing','a_iter','E_iter','bestSample','param'); 
        end
    end
end

param_cell.Cs_iter=Cs_iter;
param_cell.cc_iter=cc_iter;
param_cell.Ns_iter=Ns_iter;
param_cell.Gs_iter=Gs_iter;
param_cell.Xs_iter=Xs_iter;
param_cell.E_iter=E_iter;
param_cell.alpha_iter=a_iter;
param_cell.sigma_noise_iter=sigma_noise_iter;
param_cell.sigma_g_iter=sigma_g_iter;
param_cell.Y_true=Y_true;
param_cell.Y=Y;
param_cell.bestSample=bestSample;

end

