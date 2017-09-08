function [Y]= sample_unObs(Y, param, settings)

[D N] = size(Y);
Ytemp=sample_data(param, N, D);
Y(settings.mask_unObs==1)=Ytemp(settings.mask_unObs==1);

end