function [Lamdamudsigma]=shadow_S_Liz(Sigma,twoSigma2,mu)
cx=sqrt(2*(1-mu.^2)/pi);
Lamdamudsigma=0.5./mu.*cx.*exp(-mu.^2/twoSigma2./(1-mu.^2));

