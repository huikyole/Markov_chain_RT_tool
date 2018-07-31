% On the basis of Phase_Aerosol_MoleculeMC_BRDF_Deri_Opt.m, coded
% linearization part removed for efficient forward computation (Bob)
function [P11,P21,P33,P43,P44,omega,ksca,kext]=Phase_mat_cal(lamda0,m_r,m_i,rm,lnS,r,wtr,LR,gamma_dis,pai,tau,Ltheta,N1,Na,Nb,Nc,Nd)
% Mie scattering calculation 
% particle refractive index
%m_r=1.323;%1.44;
%m_i=9.74e-6;%0.0;%0.1;%

m_p=m_r+i*m_i;

% generation of theta (scattering angle for s1 and s2 calculation)
if gamma_dis==1 % Gamma distribution adopted in Hansen's paper
   nr=r.^[(1-3*lnS)/lnS].*exp(-r/rm/lnS);%symbol relation: "rm" here = "a" in Hansen, "lns" here = "b" in Hansen 
else           % Log-distribution used in MISR (Mike) -> Abdou 1997 paper and Grainger 2004 paper (Eq.(49))  
   kt=-[log(r)-log(rm)].^2/(2*lnS^2);
   nr=1/sqrt(2*pi)/lnS./r.*exp(kt);
end

[F11,F21,F33,F43,Qsca,Qext]=Sphere_Mie_Deri_Opt(lamda0,m_p,r,pai,tau,LR,N1,Na,Nb,Nc,Nd);        

pir2=pi*r.^2;
Interg_Core=pir2.*Qsca.*nr;
ksca=sum(Interg_Core.*wtr);

Interg_Core=pir2.*Qext.*nr;
kext=sum(Interg_Core.*wtr);     
omega=ksca/kext; %albedo Eq.(2.50) 

Fa11(Ltheta)=0;
Fa21(Ltheta)=0;
Fa33(Ltheta)=0;
Fa43(Ltheta)=0;
nrwtr=nr.*wtr;
for j=1:Ltheta
    Fa11(j)=sum(conj(F11(:,j)').*nrwtr); 
    Fa21(j)=sum(conj(F21(:,j)').*nrwtr);                
    Fa33(j)=sum(conj(F33(:,j)').*nrwtr);          
    Fa43(j)=sum(conj(F43(:,j)').*nrwtr);           
end
Fc=4*pi/(2*pi/lamda0)^2;
P11=Fc*Fa11/ksca;
P21=Fc*Fa21/ksca;
P33=Fc*Fa33/ksca;    
P43=Fc*Fa43/ksca; 
P44=P33;
return

% A copy from function "Sphere_Mie_Deri", code optimized
function [F11,F21,F33,F43,Qsca,Qext]=Sphere_Mie_Deri_Opt(lamda0,m_in,R_in,pai,tau,LR,Nmax,Na,Nb,Nc,Nd)
%tic 
m=m_in;
% Following 4 sentence of allocate space first is important for ensure computation efficiency        
anM(LR,Nmax)=0;    
bnM(LR,Nmax)=0;       
% ----------------
an(Nmax)=0;
bn(Nmax)=0;


for j=1:length(R_in)
    
    r=R_in(j);
    x=(pi*2*r/lamda0);
    y=m*x;
    % downward recursion: following n_step refer to Eq.(A11) of Ref. [3] which is comparable to Eq.(38)-part3 of Ref. [1]
    %N1=ceil(x+7.5*x^(0.34)+2);
    %N1=ceil(x+4.05*x.^(1/3)+8);
    N1=ceil(Na*x+Nb*x^Nc+Nd);
    % Lentz method to calculate Ln in Ref. [1], "Dn" in Eq.(A12) pf Ref. [3] and "An" in Eq.(36) of Ref.[2]
    Ly(N1)=fai(y,N1);
    Lx(N1)=fai(x,N1);
    n=N1;
    while n>1
        Ly(n-1)=n/y-1/[n/y+Ly(n)]; % Eq.(27) of Ref.[1]
        Lx(n-1)=n/x-1/[n/x+Lx(n)];   
        n=n-1;
    end
    % (1): n=1
    n=1;

        Fxn=2/(x^2)*(2*n+1);
        noverx=n/x;
        Lyoverm=Ly(n)/m;
        mLy=m*Ly(n);

    A(n)=1/[1-i*(cos(x)+x*sin(x))/(sin(x)-x*cos(x))];
    B0=i;
    B(n)=-noverx+1/[noverx-B0];
    Ta=[Lyoverm-Lx(n)]/[Lyoverm-B(n)];
    Tb=[mLy-Lx(n)]/[mLy-B(n)];
    an(n)=A(n)*Ta; %Eq.(28) of Ref.[1]
    bn(n)=A(n)*Tb; %Eq.(29) of Ref.[1]
    anM(j,n)=(2*n+1)/n/(n+1)*an(n);    
    bnM(j,n)=(2*n+1)/n/(n+1)*bn(n);  
    kext=Fxn*real(an(n)+bn(n));
    ksca=Fxn*(an(n)*conj(an(n))+bn(n)*conj(bn(n)));
      
    
    % (2): n>=1
    for n=2:N1
        Fxn=2/(x^2)*(2*n+1);
        nm1=n-1;
        noverx=n/x;
        Lyoverm=Ly(n)/m;
        mLy=m*Ly(n);
        B(n)=-noverx+1/[noverx-B(nm1)];
        A(n)=A(nm1)*[B(n)+noverx]/[Lx(n)+noverx];     %Eq.(34) of Ref.[1]
        Ta=[Lyoverm-Lx(n)]/[Lyoverm-B(n)];
        Tb=[mLy-Lx(n)]/[mLy-B(n)];   
        an(n)=A(n)*Ta; %Eq.(28) of Ref.[1]
        bn(n)=A(n)*Tb; %Eq.(29) of Ref.[1]
        anM(j,n)=(2*n+1)/n/(n+1)*an(n);    
        bnM(j,n)=(2*n+1)/n/(n+1)*bn(n);                  
        kext=kext+Fxn*real(an(n)+bn(n));
        ksca=ksca+Fxn*(an(n)*conj(an(n))+bn(n)*conj(bn(n)));
 
    end
    
    Qext(j)=kext;Qsca(j)=ksca;

end

s1=anM*pai(1:Nmax,:)+bnM*tau(1:Nmax,:);
s2=anM*tau(1:Nmax,:)+bnM*pai(1:Nmax,:);


F11=0.5*(s1.*conj(s1)+s2.*conj(s2));
F21=0.5*(s2.*conj(s2)-s1.*conj(s1));%Hansen Eq.(2.41) in the 1974 review article is wrong for F21 (see Liou book Eq.5.2.111b) for correct results
F33=0.5*(s1.*conj(s2)+s2.*conj(s1));
F43=0.5*i*(s1.*conj(s2)-s2.*conj(s1));

%toc
return

 
function[flk]=fai(Z,n)
%calculation of a1 when k=1 --- Eq.(15) of Shen's 1997 J. USST paper
a1=(2*n+1)/Z;
%calculation of a2 when k=2 --- Eq.(15) of Shen's 1997 J. USST paper
a2=-(2*n+3)/Z;
flk1=a1;
df=a2+1/a1;
fll=a2;
k=2;
while k
   flk=flk1*df/fll;
   if abs(flk-flk1)<1.0e-5
      break
   end
   a3=(-1)^(k+2)*[(2*(n+k+1)-1)/Z];   
   df=a3+1/df;
   fll=a3+1/fll;
   flk1=flk;
   k=k+1;
end
flk=-n/Z+flk;
return