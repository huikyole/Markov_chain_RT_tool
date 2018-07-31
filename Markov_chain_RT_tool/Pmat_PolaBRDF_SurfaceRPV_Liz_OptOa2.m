% This subroutine is copied from Pmat_PolaBRDF_SurfaceRPV_Liz_Opt
function [P11,P12,P22,P33,P34,P44]=Pmat_PolaBRDF_SurfaceRPV_Liz_OptOa2(n,v,epsirol,alamda,b,k,mu0,mu,cosfaipfai0,faiL)
% n=1.5;
% epsirol=1.0;%0.8; % epsirol is typically between [0,1]
% Sigma=0.0800;

% b=-0.5;
% k=4.5;
% alamda=0.17096;

%mu0=0.5;
%mu=0.25;

%faipfai0=3*pi/4;

cosSA=-mu*mu0+sqrt(1-mu.^2)*sqrt(1-mu0^2).*cosfaipfai0;

%1: RPV part 
f = 1/pi*(mu*mu0.*(mu+mu0)).^(k-1)*alamda.*exp(b*cosSA);% Diner's expression without hot-spot inclusion

%2: polarization part
omega=acos(cosSA); %0<omega<pi
gamma=(pi-omega)/2; % gamma>0

cosgamma=cos(gamma);
singamma=sin(gamma);
singammap=1/n*singamma;
cosgammap=sqrt(1-singammap.^2);        

cosbeta=(mu+mu0)./(2*cosgamma); %Eq.(4-5)

%(1): uniform distribution of facets  % Eq.(4-16)
%ppbeta(1:faiL)=1/2/pi; %pi is divided out of the loop
%(2): for a surface covered by hemispheres of varyinug radii % Eq.(4-17)
%ppbeta=cosbeta/pi; %pi is divided out of the loop
%(3): Gaussian facet distribution (for rough ocean surface)
%v=2; %surface wind speed (m/s)
%twoSigma2=0.003+0.00512*v;%This formulas refers to Fan et al, An improved atmospheric vector RT model incorporating rough ocean boundaries,2010
%Sigma=sqrt(twoSigma2/2);

twoSigma2=0.003+0.00512*v;%This formulas refers to Fan et al, An improved atmospheric vector RT model incorporating rough ocean boundaries,2010
Sigma=sqrt(twoSigma2/2);
tanbeta2=(1-cosbeta.^2)./(cosbeta.^2);
x=-tanbeta2/twoSigma2;
prefactor=1/twoSigma2./cosbeta.^3/pi; 
ppbeta=prefactor.*exp(x);  % Eq.(4-15)

ppbetads=-1/pi/Sigma^3./cosbeta.^3.*exp(x).*(1-tanbeta2/twoSigma2);

% to prevent that P11=0 (Pmat=0)
if ppbeta==0
    ppbeta=1e-100;
end

%Taking the shadowing effect of the sea surface wave into
%account (refers to Fan et al, An improved atmospheric vector RT model incorporating rough ocean boundaries,2010))
Lamdamu=shadow_S(Sigma,twoSigma2,mu);
Lamdamuds=shadow_S_Liz(Sigma,twoSigma2,mu);
Lamdamu0=shadow_S(Sigma,twoSigma2,mu0);           
Lamdamu0ds=shadow_S_Liz(Sigma,twoSigma2,mu0);
S=1./(1+Lamdamu+Lamdamu0); %Eq.(14) of Fan et al's paper
Sds=-(Lamdamuds+Lamdamu0ds)./(1+Lamdamu+Lamdamu0).^2;


rp=(n*cosgamma-cosgammap)./(n*cosgamma+cosgammap); % Eq.(4-7)
rs=(cosgamma-n*cosgammap)./(cosgamma+n*cosgammap);


F11=1/2*(abs(rp).^2+abs(rs).^2); % Eq. (4-8)
F12=1/2*(abs(rp).^2-abs(rs).^2);
F33=1/2*(rp.*conj(rs)+conj(rp).*rs); 
F34=i/2*(conj(rp).*rs-rp.*conj(rs));  


Fcf=ppbeta./(4*cosbeta.*mu*mu0);

%Fp1=mu;
%Fp2=mu*epsirol.*S;

% consider the "constF=pi./XMUJ/2", the above 2 sentenses are reivised as
Fp1=pi/2;
Fp2=pi/2*epsirol.*S;


    


f2=Fcf.*F11; %Eq.(4.14)
f1=f;

P11=Fp1.*f1+Fp2.*f2;

% F22=F11
P22=Fp2.*f2;


f2=Fcf.*F12; %Eq.(4.14)
P12=Fp2.*f2;

% F21=F12
%P21=P12;
   

f2=Fcf.*F33; %Eq.(4.14)
P33=Fp2.*f2;

% F33=F44
P44=P33;



f2=Fcf.*F34; %Eq.(4.14)
P34=Fp2.*f2;



end