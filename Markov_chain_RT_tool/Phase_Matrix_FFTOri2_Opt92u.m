% A combination of Phase_Matrix_FFTOri2_Opt9u.m &
% Phase_Matrix_FFTOri2_Opt2Olga2u.m for alternative choice of DStokes=3 or 4
function [PR_M,PT_M,PR_Mas,PT_Mas]=Phase_Matrix_FFTOri2_Opt92u(Lmiu0,Lmiue,M_Max,L_max,xi,Wxi,theta,P11_nr,P22_nr,P21_nr,P33_nr,P43_nr,P44_nr,P0n_bar,P2n_bar,R2n_bar,T2n_bar,Tm_bare,Rm_bare,Pnm_bare,Tm_bar0,Rm_bar0,Pnm_bar0,DStokes)
%XMU is within [-1,0] and is Gaussian intergral points in the main prog


%PhaseFun_SizeDistributionMC
%load xuas1

%M_Max=25;    
%L_max=30;
% if L_max<M_Max
%     error('Please set L_max >= M_Max !!')
% end



% B: Legendre function
% Using Liou's Eqs.(E9)&(E7) to calculate pn=Pnm(m=0)
a1=interp1(cos(theta),P11_nr,xi,'spline');
b1=interp1(cos(theta),P21_nr,xi,'spline');
a3=interp1(cos(theta),P33_nr,xi,'spline');
b2=interp1(cos(theta),-P43_nr,xi,'spline'); % Minus sign before F43_nr due to differnt sign used for P-matrix in Hansen and Siewert's papers
%a2=a1; % only for sphere    
a2=interp1(cos(theta),P22_nr,xi,'spline');
%a4=a3; % no longer valid when Rayleigh scattering is integrated into the aerosol phase function (see Eq.(2.15)-(2.16) in Hansen 1974's and the subroutine Phase_Aerosol_Molecule.m)
a4=interp1(cos(theta),P44_nr,xi,'spline');



F2L1=2*[0:1:L_max]+1;
ones_arr=ones(Lmiue,1);

% here to be added for Olga (2)
Lxi=length(xi);
alfa1(L_max+1)=0;
alfa2(L_max+1)=0;
alfa3(L_max+1)=0;
alfa4(L_max+1)=0;
beta1(L_max+1)=0;
beta2(L_max+1)=0;
R2n_bar(L_max+2,Lxi)=0;
T2n_bar(L_max+2,Lxi)=0;
% =============================

% Part 1: calculation of basic constants (alfa, beta, gamma, zeta, delta, epsirol)
% a) beta & delta
for L=0:L_max    
    alfa1(L+1)=(L+0.5)*sum(Wxi.*a1.*P0n_bar(L+1,:));
    alfa4(L+1)=(L+0.5)*sum(Wxi.*a4.*P0n_bar(L+1,:));
end

% b) gamma & epsirol
for L=2:L_max    
    beta1(L+1)=(L+0.5)*sum(Wxi.*b1.*P2n_bar(L+1,:)); 
    beta2(L+1)=(L+0.5)*sum(Wxi.*b2.*P2n_bar(L+1,:));
end



% c) zeta & alfa
m=2;  
for L=m:L_max % calculate for L
    alfa3(L+1)=(L+0.5)*sum(Wxi.*(a3.*R2n_bar(L+1,:)+a2.*T2n_bar(L+1,:)));
    alfa2(L+1)=(L+0.5)*sum(Wxi.*(a2.*R2n_bar(L+1,:)+a3.*T2n_bar(L+1,:)));
end


alfa1=repmat(alfa1,Lmiue,1);
alfa2=repmat(alfa2,Lmiue,1);
alfa3=repmat(alfa3,Lmiue,1);
alfa4=repmat(alfa4,Lmiue,1);
beta1=repmat(beta1,Lmiue,1);
beta2=repmat(beta2,Lmiue,1);

PR_M(DStokes*Lmiue/2,DStokes*Lmiu0/2,M_Max+1)=0;
PT_M(DStokes*Lmiue/2,DStokes*Lmiu0/2,M_Max+1)=0;
PR_Mas(DStokes*Lmiue/2,DStokes*Lmiu0/2,M_Max+1)=0;
PT_Mas(DStokes*Lmiue/2,DStokes*Lmiu0/2,M_Max+1)=0;


Lmiued2_arr=[1:Lmiue/2];
Lmiu0d2_arr=[1:Lmiu0/2];

for m=0:M_Max
    Am11=Pnm_bare(m+1:L_max+1,:,m+1)'.*alfa1(:,m+1:L_max+1)*Pnm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
    share1=Pnm_bare(m+1:L_max+1,:,m+1)'.*beta1(:,m+1:L_max+1);
    Am12=share1*Rm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
    Am13=-share1*Tm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
    %Am14=0;

    Am21=Rm_bare(m+1:L_max+1,:,m+1)'.*beta1(:,m+1:L_max+1)*Pnm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
    share1=Rm_bare(m+1:L_max+1,:,m+1)'.*alfa2(:,m+1:L_max+1);
    share2=Tm_bare(m+1:L_max+1,:,m+1)'.*alfa3(:,m+1:L_max+1);
    Am22=share1*Rm_bar0(m+1:L_max+1,1:Lmiu0,m+1)+share2*Tm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
    Am23=-(share1*Tm_bar0(m+1:L_max+1,1:Lmiu0,m+1)+share2*Rm_bar0(m+1:L_max+1,1:Lmiu0,m+1));

    Am31=-Tm_bare(m+1:L_max+1,:,m+1)'.*beta1(:,m+1:L_max+1)*Pnm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
    share1=Tm_bare(m+1:L_max+1,:,m+1)'.*alfa2(:,m+1:L_max+1);
    share2=Rm_bare(m+1:L_max+1,:,m+1)'.*alfa3(:,m+1:L_max+1);
    Am32=-(share1*Rm_bar0(m+1:L_max+1,1:Lmiu0,m+1)+share2*Tm_bar0(m+1:L_max+1,1:Lmiu0,m+1));
    Am33=share1*Tm_bar0(m+1:L_max+1,1:Lmiu0,m+1)+share2*Rm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
  
    
    mm=1;nn=1;
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am11(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am11(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am11(Lmiue/2+1:Lmiue,1:Lmiu0/2));
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am11(1:Lmiue/2,1:Lmiu0/2)));    

    mm=1;nn=2;
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am12(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am12(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);   
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am12(Lmiue/2+1:Lmiue,1:Lmiu0/2));    
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am12(1:Lmiue/2,1:Lmiu0/2)));    
    
    mm=1;nn=3;      
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-flipud(Am13(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-Am13(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-fliplr(Am13(Lmiue/2+1:Lmiue,1:Lmiu0/2));    
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-flipud(fliplr(Am13(1:Lmiue/2,1:Lmiu0/2)));
    
    
    mm=2;nn=1;    
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am21(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am21(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);    
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am21(Lmiue/2+1:Lmiue,1:Lmiu0/2));
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am21(1:Lmiue/2,1:Lmiu0/2)));    
    
    mm=2;nn=2;      
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am22(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am22(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am22(Lmiue/2+1:Lmiue,1:Lmiu0/2));
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am22(1:Lmiue/2,1:Lmiu0/2)));    
    
    mm=2;nn=3;      
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-flipud(Am23(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-Am23(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);    
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-fliplr(Am23(Lmiue/2+1:Lmiue,1:Lmiu0/2));    
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-flipud(fliplr(Am23(1:Lmiue/2,1:Lmiu0/2)));

    
    mm=3;nn=1;     
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am31(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am31(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am31(Lmiue/2+1:Lmiue,1:Lmiu0/2));    
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am31(1:Lmiue/2,1:Lmiu0/2)));    
    
    mm=3;nn=2;      
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am32(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am32(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);    
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am32(Lmiue/2+1:Lmiue,1:Lmiu0/2));
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am32(1:Lmiue/2,1:Lmiu0/2)));

    mm=3;nn=3;      
    PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am33(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
    PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am33(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);    
    PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am33(Lmiue/2+1:Lmiue,1:Lmiu0/2));
    PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am33(1:Lmiue/2,1:Lmiu0/2)));      

end


if DStokes==4
    for m=0:M_Max
        Am24=-Tm_bare(m+1:L_max+1,:,m+1)'.*beta2(:,m+1:L_max+1)*Pnm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
        Am34=Rm_bare(m+1:L_max+1,:,m+1)'.*beta2(:,m+1:L_max+1)*Pnm_bar0(m+1:L_max+1,1:Lmiu0,m+1); 
        
        share1=Pnm_bare(m+1:L_max+1,:,m+1)'.*beta2(:,m+1:L_max+1);
        Am42=share1*Tm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
        Am43=-share1*Rm_bar0(m+1:L_max+1,1:Lmiu0,m+1);
        Am44=Pnm_bare(m+1:L_max+1,:,m+1)'.*alfa4(:,m+1:L_max+1)*Pnm_bar0(m+1:L_max+1,1:Lmiu0,m+1);                  
        
        mm=2;nn=4;     
        PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-flipud(Am24(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
        PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-Am24(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);    
        PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-fliplr(Am24(Lmiue/2+1:Lmiue,1:Lmiu0/2));
        PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=-flipud(fliplr(Am24(1:Lmiue/2,1:Lmiu0/2)));  

        mm=3;nn=4;     
        PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am34(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
        PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am34(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);
        PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am34(Lmiue/2+1:Lmiue,1:Lmiu0/2));
        PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am34(1:Lmiue/2,1:Lmiu0/2)));   

        mm=4;nn=2;     
        PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am42(1:Lmiue/2,1+Lmiu0/2:Lmiu0));             
        PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am42(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);    
        PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am42(Lmiue/2+1:Lmiue,1:Lmiu0/2));    
        PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am42(1:Lmiue/2,1:Lmiu0/2)));

        mm=4;nn=3;     
        PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am43(1:Lmiue/2,1+Lmiu0/2:Lmiu0));
        PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am43(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0);
        PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am43(Lmiue/2+1:Lmiue,1:Lmiu0/2));
        PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am43(1:Lmiue/2,1:Lmiu0/2)));

        mm=4;nn=4;     
        PR_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(Am44(1:Lmiue/2,1+Lmiu0/2:Lmiu0));            
        PT_M(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=Am44(1+Lmiue/2:Lmiue,1+Lmiu0/2:Lmiu0); 
        PR_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=fliplr(Am44(Lmiue/2+1:Lmiue,1:Lmiu0/2));     
        PT_Mas(DStokes*(Lmiued2_arr-1)+mm,DStokes*(Lmiu0d2_arr-1)+nn,m+1)=flipud(fliplr(Am44(1:Lmiue/2,1:Lmiu0/2)));   
    end
end
    
return

