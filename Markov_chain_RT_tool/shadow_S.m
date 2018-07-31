function [Lamdamu]=shadow_S(Sigma,twoSigma2,mu)
t1=sqrt(2*(1-mu.^2)/pi);
t2=Sigma./mu.*exp(-mu.^2/twoSigma2./(1-mu.^2));
t3=erfc(mu./Sigma./sqrt(2*(1-mu.^2)));
Lamdamu=0.5*(t1.*t2-t3);
