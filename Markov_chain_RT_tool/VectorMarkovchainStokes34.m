% Based on HomogeneousAtmosphereV_O~12, constructed as a compelte package of Forward model 
function [IS,IM]=VectorMarkovchainStokes34(deta,b1,k1,r01,nn,vv,epsirol,M_Max,XMU,WTMU,XMU0ALL,XMUeALL,faipfai0,Lfaipfai0,NUMLAY_GROUP,NMultiplication,L_max,xi_number,DStokes)
global Aerosol_Type SSA_p_TYPE MiemixtypeforL Miemixtype LMiemixtype AF_p theta TAUALL RAYF OMEGA P11_nm_TYPE P21_nm_TYPE P22_nm_TYPE P33_nm_TYPE P43_nm_TYPE P44_nm_TYPE
tic  
DStokes1=DStokes-1;
LXMU0ALL=length(XMU0ALL);
LXMUeALL=length(XMUeALL);

%Following 3 sentences: the last layer does not have effect, it is actually for surface
TAUALL=[TAUALL,TAUALL(end)*1.5]; 
RAYF=[RAYF,1]; 
MIEF=1-RAYF;
OMEGA=[OMEGA,1];


% sub-group division
NUMGROUP=floor(length(RAYF)/NUMLAY_GROUP);
GROUPLAYER=[1:NUMLAY_GROUP:NUMLAY_GROUP*NUMGROUP];
if mod(length(RAYF),NUMLAY_GROUP)~=0
    GROUPLAYER=[GROUPLAYER,length(RAYF)]; %Note that here orignial "42" changes to "43" which indicates the last layer having no effect but for surface 
end
% number for performing "adding-algorithm"
NUMADD=length(GROUPLAYER)-1;

% Rayleigh scattering   
% depolarization factor: 0 for isotropic Rayleigh scattering
%deta=0.0; %  2.908633816103590d-2;% 
DDeta=(1-deta)/(1+deta/2);
DDetap=(1-2*deta)/(1-deta);  
P11_nr=DDeta*3/4*(1+cos(theta).^2)+(1-DDeta)*1;
P22_nr=DDeta*3/4*(1+cos(theta).^2);
P21_nr=-DDeta*3/4*sin(theta).^2+0;
P33_nr=DDeta*3/2*cos(theta)+0;
P43_nr(1:length(theta))=0;
argue_val=0; %0 in Hansen's paper (also Diner's opinion) while 1 in Emde's paper 
P44_nr=DDeta*DDetap*3/2*cos(theta)+(1-DDeta)*argue_val;     

%% Following for checking the scalar code
% P22_nr(1:end)=0;
% P21_nr(1:end)=0;
% P33_nr(1:end)=0;
% P43_nr(1:end)=0;
% P44_nr(1:end)=0;

% Inci_Ang_Location=[length(XMU)+1]/2; %check to see whether this is equal to 60 deg.
% if abs(acos(XMU(Inci_Ang_Location))*180/pi-60)>1.0e-10
%     display('The incidence agle 60 deg is not within the Gaussian points! ')
%     return
% end

%tic
% FFT series of the depolarizing BRDF surface (since only P11 exists in the phase matrix, the FFT series can evaluated in the same way as done for scalar case)
LMU=length(XMU);
XMV=sqrt(1-XMU.^2);
PHI_number=1500;
%[PHI,WTPHI]=Angle_Gauss_10pointMC(0,pi,PHI_number); 
[PHI,WTPHI]=lgwt(PHI_number,0,2*pi);
WTPHI=WTPHI/pi; % for consistence with Espoito 
WTPHI=repmat(WTPHI,LXMUeALL,1);
m_arr=[0:1:M_Max]';
mPHI=(m_arr*PHI)';
 
%allocate the space first to ensure the computation efficiency
PRSur_M(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1)=0;
if (r01~=0)|(epsirol~=0)
    for ni=1:LXMU0ALL
        %for nj=1:LXMUeALL
            XMUI=XMU0ALL(ni);XMUJ=XMUeALL;%(nj);
            XMVI=sqrt(1-XMUI^2);XMVJ=sqrt(1-XMUJ.^2);
            XIXJ=XMUI*XMUJ;VIVJ=XMVI*XMVJ;
            CSTHR=-repmat(XIXJ',1,PHI_number)+VIVJ'*cos(PHI);

            XMUJmat=repmat(XMUJ',1,PHI_number);
            XMVJmat=repmat(XMVJ',1,PHI_number); 
            PHImat=repmat(PHI,LXMUeALL,1);
            cosi1R=[XMUJmat+XMUI*CSTHR]./sqrt(1-CSTHR.^2)/XMVI;        
            cosi2R=[-XMUI-XMUJmat.*CSTHR]./sqrt(1-CSTHR.^2)./XMVJmat; 
            sini1R=[XMVJmat.*sin(-PHImat)]./sqrt(1-CSTHR.^2); 
            sini2R=[XMVI*sin(-PHImat)]./sqrt(1-CSTHR.^2);             
            cos2alfa0=2*cosi1R.^2-1;
            sin2alfa0=2*sini1R.*cosi1R;

            cos2alfa=2*cosi2R.^2-1;
            sin2alfa=2*sini2R.*cosi2R;  



            % without polarization with polarization
            cosPHI=cos(PHImat); %use "-PHI" to get consistence with WENG's Eq.(20) where "fai0-fai" is used
            sinPHI=sin(PHImat); %use "-PHI" to get consistence with WENG's Eq.(20) where "fai0-fai" is used 

            [P11,P12,P22,P33,P34,P44]=Pmat_PolaBRDF_SurfaceRPV_Liz_OptOa2(nn,vv,epsirol,r01,b1,k1,XMUI,XMUJmat,cosPHI,PHI_number);


            [RPRsur11,RPRsur12,RPRsur13,...
             RPRsur21,RPRsur22,RPRsur23,RPRsur24,...
             RPRsur31,RPRsur32,RPRsur33,RPRsur34,...
             RPRsur42,RPRsur43,RPRsur44]=RPR_Cal(P11,P12,P22,P33,P34,P44,cos2alfa0,sin2alfa0,cos2alfa,sin2alfa);

            %constF=pi./XMUJ/2; % this sentense is put into the subroutine "Pmat_PolaBRDF_SurfaceRPV_Liz_OptOa2"


            % merge the consine and sine part together to save time in calculation
            PRSur_M(1:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+1,1:M_Max+1)=RPRsur11.*WTPHI*[cos(mPHI)+sin(mPHI)];
            PRSur_M(1:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+2,1:M_Max+1)=RPRsur12.*WTPHI*[cos(mPHI)+sin(mPHI)];
            PRSur_M(1:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+3,1:M_Max+1)=RPRsur13.*WTPHI*[cos(mPHI)+sin(mPHI)];
            PRSur_M(2:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+1,1:M_Max+1)=RPRsur21.*WTPHI*[cos(mPHI)+sin(mPHI)];
            PRSur_M(2:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+2,1:M_Max+1)=RPRsur22.*WTPHI*[cos(mPHI)+sin(mPHI)];
            PRSur_M(2:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+3,1:M_Max+1)=RPRsur23.*WTPHI*[cos(mPHI)+sin(mPHI)];
            PRSur_M(3:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+1,1:M_Max+1)=RPRsur31.*WTPHI*[cos(mPHI)+sin(mPHI)];
            PRSur_M(3:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+2,1:M_Max+1)=RPRsur32.*WTPHI*[cos(mPHI)+sin(mPHI)];
            PRSur_M(3:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+3,1:M_Max+1)=RPRsur33.*WTPHI*[cos(mPHI)+sin(mPHI)];

            if DStokes==4
                PRSur_M(4:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+2,1:M_Max+1)=RPRsur42.*WTPHI*[cos(mPHI)+sin(mPHI)];
                PRSur_M(4:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+3,1:M_Max+1)=RPRsur43.*WTPHI*[cos(mPHI)+sin(mPHI)];

                PRSur_M(2:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+4,1:M_Max+1)=RPRsur24.*WTPHI*[cos(mPHI)+sin(mPHI)];
                PRSur_M(3:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+4,1:M_Max+1)=RPRsur34.*WTPHI*[cos(mPHI)+sin(mPHI)];            
                PRSur_M(4:DStokes:LXMUeALL*DStokes,(ni-1)*DStokes+4,1:M_Max+1)=RPRsur44.*WTPHI*[cos(mPHI)+sin(mPHI)];        
            end

        %end
    end  
end

%toc

%L_max=1313; % Moment number for FFT decomposition of phase matrix
if L_max<M_Max
    L_max=M_Max; % since L_max should be >=M_Max
end
%convergence of the numerical result of alfa, beta, zeta and epsirol.
% xi_number=1500;  % check to see whether this value is large enough to ensure the convergence of the numerical integral.
                % should be as large as 30 for 1um particle 60 for 9um
                % particle size
xi_min=-1; %check
xi_max=1;
%xi_number=30;%15; %should be as large as 60 for 9um particle size 
Lxi=xi_number;
[xi,Wxi]=lgwt(Lxi,xi_min,xi_max);

% here to be added for Oa (1)
P0n(1:L_max+1,1:Lxi)=0;
P2n(1:L_max+1,1:Lxi)=0;
P0n_bar(1:L_max+1,1:Lxi)=0;
P2n_bar(1:L_max+1,1:Lxi)=0;
R2n_bar(L_max+2,Lxi)=0;
T2n_bar(L_max+2,Lxi)=0;
% =============================

% for m=0;
P0n(1,1:Lxi)=1; %n=0
P0n(2,1:Lxi)=xi; %n=1
for n=2:L_max
    P0n(n+1,1:Lxi)=[(2*(n-1)+1)*xi.*P0n(n,1:Lxi)-(n-1)*P0n(n-1,1:Lxi)]/(n-1+1); %Eq.(E7) of Liou
end 
% normalizaition by [(n-m)!(n+m)!]*(1/2)
m=0;
for n=m:L_max
    P0n_bar(n+1,:)=P0n(n+1,:)/sqrt(prodc(n-m,n+m));
end
% for m=2;
P2n(3,1:Lxi)=3*(1-xi.^2); %n=2
P2n(4,1:Lxi)=15*xi.*(1-xi.^2); %n=3
for n=4:L_max
    P2n(n+1,1:Lxi)=[(2*(n-1)+1)*xi.*P2n(n,1:Lxi)-(n-1+2)*P2n(n-1,1:Lxi)]/(n-1-1); %Eq.(21) of Siewert
end
% normalizaition by [(n-m)!(n+m)!]*(1/2)
m=2;
for n=m:L_max
    P2n_bar(n+1,:)=P2n(n+1,:)/sqrt(prodc(n-m,n+m));
end

m=2;  
L=m;
R2n_bar(L+1,:)=sqrt(m*(m-1)/(m+1)/(m+2))*(1+xi.^2)./(1-xi.^2).*P2n_bar(L+1,:);  
T2n_bar(L+1,:)=sqrt(m*(m-1)/(m+1)/(m+2))*(2*xi)./(1-xi.^2).*P2n_bar(L+1,:);
for L=m:L_max % then calculate for L+1
    tmL=sqrt(((L+1)^2-m^2)*((L+1)^2-4))/(L+1);
    if L==m
        R2n_bar(L+1+1,:)=(2*L+1)*(xi.*R2n_bar(L+1,:)-2*m/L/(L+1)*T2n_bar(L+1,:))/tmL;  
        T2n_bar(L+1+1,:)=(2*L+1)*(xi.*T2n_bar(L+1,:)-2*m/L/(L+1)*R2n_bar(L+1,:))/tmL;         
    else
        R2n_bar(L+1+1,:)=(2*L+1)*(xi.*R2n_bar(L+1,:)-2*m/L/(L+1)*T2n_bar(L+1,:))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*R2n_bar(L-1+1,:)/tmL;  
        T2n_bar(L+1+1,:)=(2*L+1)*(xi.*T2n_bar(L+1,:)-2*m/L/(L+1)*R2n_bar(L+1,:))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*T2n_bar(L-1+1,:)/tmL;   
    end
end
                 
% 1) for the scattering angle
miue=[-fliplr(XMUeALL),XMUeALL];
[Tm_bare,Rm_bare,Pnm_bare]=TR(miue,L_max,M_Max);
% 2) for the incident angle
miu0=[-fliplr(XMU0ALL),XMU0ALL];
[Tm_bar0,Rm_bar0,Pnm_bar0]=TR(miu0,L_max,M_Max);
Lmiu0=2*LXMU0ALL;
Lmiue=2*LXMUeALL;
   

[PR_MRAY,PT_MRAY,PR_MasRAY,PT_MasRAY]=Phase_Matrix_FFTOri2_Opt92u(Lmiu0,Lmiue,3,L_max,xi,Wxi,theta,P11_nr,P22_nr,P21_nr,P33_nr,P43_nr,P44_nr,P0n_bar,P2n_bar,R2n_bar,T2n_bar,Tm_bare,Rm_bare,Pnm_bare,Tm_bar0,Rm_bar0,Pnm_bar0,DStokes);
PR_MRAY(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1)=0;
PT_MRAY(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1)=0;
PR_MasRAY(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1)=0;
PT_MasRAY(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1)=0;
    
Ltheta=length(theta);
P11_nm(LMiemixtype,Ltheta)=0; 
P21_nm(LMiemixtype,Ltheta)=0;
P22_nm(LMiemixtype,Ltheta)=0;    
P33_nm(LMiemixtype,Ltheta)=0;
P43_nm(LMiemixtype,Ltheta)=0;
P44_nm(LMiemixtype,Ltheta)=0; 

PR_MMIE(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1,LMiemixtype)=0;
PT_MMIE(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1,LMiemixtype)=0;
PR_MasMIE(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1,LMiemixtype)=0;
PT_MasMIE(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1,LMiemixtype)=0;

PR_MMIE_TYPE(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1,Aerosol_Type)=0;
PT_MMIE_TYPE(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1,Aerosol_Type)=0;
PR_MasMIE_TYPE(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1,Aerosol_Type)=0;
PT_MasMIE_TYPE(LXMUeALL*DStokes,LXMU0ALL*DStokes,M_Max+1,Aerosol_Type)=0;


MiemixtypeforL(end+1)=MiemixtypeforL(end); % "+1" for allocate space for surface layer

% for Mie aerosols
%% Phase matrix decomposition for all types of aerosols

for n=1:Aerosol_Type
    if SSA_p_TYPE(n)~=0
        if DStokes==4
            [PR_MMIE_TYPE(:,:,:,n),PT_MMIE_TYPE(:,:,:,n),PR_MasMIE_TYPE(:,:,:,n),PT_MasMIE_TYPE(:,:,:,n)]=Phase_Matrix_FFTOri2_Opt92u(Lmiu0,Lmiue,M_Max,L_max,xi,Wxi,theta,P11_nm_TYPE(n,:),P22_nm_TYPE(n,:),P21_nm_TYPE(n,:),P33_nm_TYPE(n,:),P43_nm_TYPE(n,:),P44_nm_TYPE(n,:),P0n_bar,P2n_bar,R2n_bar,T2n_bar,Tm_bare,Rm_bare,Pnm_bare,Tm_bar0,Rm_bar0,Pnm_bar0,DStokes);       
        else
            [PR_MMIE_TYPE(:,:,:,n),PT_MMIE_TYPE(:,:,:,n),PR_MasMIE_TYPE(:,:,:,n),PT_MasMIE_TYPE(:,:,:,n)]=Phase_Matrix_FFTOri2_Opt92u(Lmiu0,Lmiue,M_Max,L_max,xi,Wxi,theta,P11_nm_TYPE(n,:),P22_nm_TYPE(n,:),P21_nm_TYPE(n,:),P33_nm_TYPE(n,:),P43_nm_TYPE(n,:),P44_nm_TYPE(n,:),P0n_bar,P2n_bar,R2n_bar,T2n_bar,Tm_bare,Rm_bare,Pnm_bare,Tm_bar0,Rm_bar0,Pnm_bar0,DStokes);
        end
    end
end

% Determine the phase matrix and its FFT decomposition in different layers
for L=1:LMiemixtype
    layernum_here=Miemixtype(L);
    for n=1:Aerosol_Type
        PR_MMIE(:,:,:,L)=PR_MMIE(:,:,:,L)+AF_p(layernum_here,n)*PR_MMIE_TYPE(:,:,:,n);
        PT_MMIE(:,:,:,L)=PT_MMIE(:,:,:,L)+AF_p(layernum_here,n)*PT_MMIE_TYPE(:,:,:,n);
        PR_MasMIE(:,:,:,L)=PR_MasMIE(:,:,:,L)+AF_p(layernum_here,n)*PR_MasMIE_TYPE(:,:,:,n);
        PT_MasMIE(:,:,:,L)=PT_MasMIE(:,:,:,L)+AF_p(layernum_here,n)*PT_MasMIE_TYPE(:,:,:,n);

        P11_nm(L,:)=P11_nm(L,:)+AF_p(layernum_here,n)*P11_nm_TYPE(n,:);
        P21_nm(L,:)=P21_nm(L,:)+AF_p(layernum_here,n)*P21_nm_TYPE(n,:);
        P22_nm(L,:)=P22_nm(L,:)+AF_p(layernum_here,n)*P22_nm_TYPE(n,:);                
        P33_nm(L,:)=P33_nm(L,:)+AF_p(layernum_here,n)*P33_nm_TYPE(n,:); 
        P43_nm(L,:)=P43_nm(L,:)+AF_p(layernum_here,n)*P43_nm_TYPE(n,:);
        P44_nm(L,:)=P44_nm(L,:)+AF_p(layernum_here,n)*P44_nm_TYPE(n,:);                 
    end
end
             
    
IM(1:DStokes*LXMUeALL,1:LXMU0ALL,1:Lfaipfai0)=0;
Im(1:DStokes*LXMUeALL,1:LXMU0ALL,1:Lfaipfai0,1:M_Max+1)=0;
Imc(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL,1:M_Max+1)=0;
Ims(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL,1:M_Max+1)=0;
WMU12=repmat(WTMU,DStokes,1);
WMU12=reshape(WMU12,DStokes*LMU,1);
XMU12=repmat(XMU,DStokes,1);
XMU12=reshape(XMU12,DStokes*LMU,1);

XMJ12mat=repmat(XMU12',DStokes*LXMUeALL,1);
XMEALL=repmat(XMUeALL,DStokes,1);
XMEALL=reshape(XMEALL,DStokes*LXMUeALL,1); 
XME12mat=repmat(XMEALL,1,DStokes*LMU);

% change 4
XMU012ALL=repmat(XMU0ALL,DStokes,1);
XMU012ALL=reshape(XMU012ALL,DStokes*LXMU0ALL,1);  
XMU0ALLmat=repmat(XMU012ALL',DStokes*LMU,1);    
XMUIALLmat=repmat(XMU12,1,DStokes*LXMU0ALL);
WMUIALLmat=repmat(WMU12,1,DStokes*LXMU0ALL);  
% the end of change 4

%change 38 (following 3 lines added)
XMUeOUT=repmat(XMEALL,1,DStokes*LXMU0ALL);
XMUiOUT=XMUeOUT;
XMU0OUT=repmat(XMU012ALL',DStokes*LXMUeALL,1);   

ones1=[];ones2=[];
ones1(1:DStokes*LXMUeALL)=0;
ones2(1:DStokes*LXMUeALL)=0;
ones1(1,3:DStokes:end)=1;
ones2(1,1:DStokes:end)=1;
ones2(1,2:DStokes:end)=1;

if DStokes==4
    ones1(1,4:DStokes:end)=1;
    onescosmatRS=repmat([ones2;ones2;ones1;ones1]',1,LXMU0ALL);
    onessinmatRS=repmat([ones1;ones1;ones2;ones2]',1,LXMU0ALL);
else
    onescosmatRS=repmat([ones2;ones2;ones1]',1,LXMU0ALL);
    onessinmatRS=repmat([ones1;ones1;ones2]',1,LXMU0ALL);      
end



SLAST(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL,1:M_Max+1)=0;
    
%tic

NLAST=GROUPLAYER(end);
for N=1:NUMADD
    N
    
    NFIRST=GROUPLAYER(NUMADD-N+1);
    NBIG=NLAST-NFIRST+1;
    if NFIRST>1
        tau=[0,TAUALL(NFIRST:NFIRST+NBIG-1)-TAUALL(NFIRST-1)];
    else
        tau=[0,TAUALL(NFIRST:NFIRST+NBIG-1)];
    end
    omega=OMEGA(NFIRST:NFIRST+NBIG-1);
    Fsr=RAYF(NFIRST:NFIRST+NBIG-1);  
    Fsm=MIEF(NFIRST:NFIRST+NBIG-1); 

   
    % for later cosine and sine components extraction from the Q-matrix
    ones1=[];ones2=[];
    ones1(1:2*DStokes*NBIG*LMU)=0;
    ones2(1:2*DStokes*NBIG*LMU)=0;
    ones1(1,3:DStokes:end)=1;
    ones2(1,1:DStokes:end)=1;
    ones2(1,2:DStokes:end)=1;

  
    if DStokes==4
        ones1(1,4:DStokes:end)=1;
        
        onescos=[ones2;ones2;ones1;ones1];
        onessin=[ones1;ones1;ones2;ones2];
        onescosmatS=repmat([ones2',ones2',ones1',ones1'],1,LXMU0ALL);
        onessinmatS=repmat([ones1',ones1',ones2',ones2'],1,LXMU0ALL); 
    else
        onescos=[ones2;ones2;ones1];
        onessin=[ones1;ones1;ones2];
        onescosmatS=repmat([ones2',ones2',ones1'],1,LXMU0ALL);% change 31
        onessinmatS=repmat([ones1',ones1',ones2'],1,LXMU0ALL);% change 31         
    end
    onessinmatQ=repmat(onessin,2*NBIG*LMU,1);
    onescosmatR=repmat(onescos,LXMUeALL,1);% change 30 
    onessinmatR=repmat(onessin,LXMUeALL,1);% change 31  
    
    Qrown=2*DStokes*LMU*NBIG;
    Qrown2=Qrown/2;
    Qmkl(1:Qrown,1:Qrown,1:M_Max+1)=0;
    Rmle(1:DStokes*LXMUeALL,1:Qrown,1:M_Max+1)=0;
    Rle(1:DStokes*LXMUeALL,1:Qrown,1:M_Max+1)=0;
    Rles(1:DStokes*LXMUeALL,1:Qrown,1:M_Max+1)=0;    
    SOURCE1(1:Qrown2,1:DStokes*LXMU0ALL,1:M_Max+1)=0;
    SOURCE2(1:Qrown2,1:DStokes*LXMU0ALL,1:M_Max+1)=0;
    SOURCEs(1:Qrown,1:DStokes*LXMU0ALL,1:M_Max+1)=0;
    QPI(1:Qrown,1:DStokes*LXMU0ALL,1:M_Max+1)=0;

    
    % Part 1: calculate the transisition matrix Qkl: k-->(i,n) and l-->(j,n')
    % in calculating Qkl-matrix, note "note that ui>0" at P1061   
    

        % calculation of W_ninp
    for n=1:NBIG    
        dtn=tau(n+1)-tau(n);
        for ii=1:LMU %incidence angle        
            for np=1:n %NBIG
                dtnp=tau(np+1)-tau(np);
                %% calculation of W(n,i,n')   
                nui=(n-1)*DStokes*LMU+DStokes*(ii-1)+1;

                %(1)n>=np (upwelling) or n=np (partly recaptured when upwelling)
                FA=exp(-dtn/XMU(ii));
                tmpa(nui:nui+DStokes1,np)=XMU(ii)/dtn*(1-FA);
                FB=exp(-(tau(n)-tau(np+1))/XMU(ii));
                tmpb(nui:nui+DStokes1,np)=FB;
                if n==np
                    FC=exp(-dtn/XMU(ii));
                    tmpc(nui:nui+DStokes1,np)=1-FC;
                    W_ninp(nui:nui+DStokes1,np)=1-XMU(ii)/dtn*tmpc(nui:nui+DStokes1,np); %STAY                    


                else 
                    FC=exp(-dtnp/XMU(ii));
                    tmpc(nui:nui+DStokes1,np)=1-FC; %ESC
                    W_ninp(nui:nui+DStokes1,np)=tmpa(nui:nui+DStokes1,np).*tmpb(nui:nui+DStokes1,np).*tmpc(nui:nui+DStokes1,np); 


                end 

                %(2) switch n<-->np so that np>=n (downwelling) or n=np (partly recaptured when downwelling)                    
                [n,np]=switchf(n,np);
                [dtn,dtnp]=switchf(dtn,dtnp);

                nui=(2*NBIG-n)*DStokes*LMU+DStokes*(ii-1)+1;

                FA=exp(-dtn/XMU(ii));
                tmpa(nui:nui+DStokes1,np)=XMU(ii)/dtn*[1-FA]; 
                FB=exp(-(tau(np)-tau(n+1))/XMU(ii));
                tmpb(nui:nui+DStokes1,np)=FB; 
                if n==np
                    FC=exp(-dtn/XMU(ii));
                    tmpc(nui:nui+DStokes1,np)=1-FC;
                    W_ninp(nui:nui+DStokes1,np)=1-XMU(ii)/dtn*tmpc(nui:nui+DStokes1,np);                     

                else 
                    FC=exp(-dtnp/XMU(ii));
                    tmpc(nui:nui+DStokes1,np)=1-FC;
                    W_ninp(nui:nui+DStokes1,np)=tmpa(nui:nui+DStokes1,np).*tmpb(nui:nui+DStokes1,np).*tmpc(nui:nui+DStokes1,np);

                end                


                % switch back after job (2) is done
                [n,np]=switchf(n,np);
                [dtn,dtnp]=switchf(dtn,dtnp);   
                % job (2) done
            end % of np            
        end % of ii
    end %of n
    %W_ninp       

    % point 1 correction for adding method integration
    %if NFIRST~=GROUPLAYER(NUMADD)  %correction if it is not the first adding
    for ii=1:LMU
        nui=(NBIG-1)*DStokes*LMU+DStokes*(ii-1)+1;
        % No photon transition from surface layer to surface layer
        W_ninp(nui:nui+DStokes1,NBIG)=0;


        nui2=(2*NBIG-NBIG)*DStokes*LMU+DStokes*(ii-1)+1;
        W_ninp(nui2:nui2+DStokes1,NBIG)=0;

        % transistion from surface layer to other layer &
        %                  other layer to surface layer
        for np=1:NBIG-1 %NLAST                    
            tmpbc=exp(-(tau(NBIG)-tau(np+1))/XMU(ii));  %PNTONP                                      
            % upwelling
            FC=1-tmpc(nui:nui+DStokes1,np);
            W_ninp(nui:nui+DStokes1,np)=tmpbc*tmpc(nui:nui+DStokes1,np);


            % downwelling
            nui3=(2*NBIG-np)*DStokes*LMU+DStokes*(ii-1)+1;
            FC=1-tmpa(nui3:nui3+DStokes1,np)*(tau(np+1)-tau(np))/XMU(ii);
            W_ninp(nui3:nui3+DStokes1,NBIG)=tmpbc*tmpa(nui3:nui3+DStokes1,np); 

        end
    end                
    %end        




    % caulation of Qmkl
    for n=1:NBIG      
        for np=1:n                     
            commonterm1=0.5*omega(np)*(WMU12*W_ninp((n-1)*DStokes*LMU+1:n*DStokes*LMU,np)'); 
            commonterm2=0.5*omega(n)*(WMU12*W_ninp((2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,n)'); 

            %% calculatin of Qm_kl,k--(i,n) and l-->(j,n') 
            %% (1) upwelling to upwelling using PT
            commonterm=commonterm1;  
            for m=0:M_Max            
                Pmat=[Fsr(np)*PT_MasRAY(1:DStokes*LMU,1:DStokes*LMU,m+1)+Fsm(np)*PT_MasMIE(1:DStokes*LMU,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+np-1))];  
                Qmkl(DStokes*(np-1)*LMU+1:DStokes*np*LMU,DStokes*(n-1)*LMU+1:DStokes*n*LMU,m+1)=commonterm.*Pmat;                  

                %% (2) upwelling to downwelling using PR
                Pmat=[Fsr(np)*PR_MasRAY(1:DStokes*LMU,1:DStokes*LMU,m+1)+Fsm(np)*PR_MasMIE(1:DStokes*LMU,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+np-1))];   
                Qmkl(DStokes*(2*NBIG-np)*LMU+1:DStokes*(2*NBIG-np+1)*LMU,DStokes*(n-1)*LMU+1:DStokes*n*LMU,m+1)=commonterm.*Pmat;                  
            end


            %(2) switch n<-->np so that np>=n (downwelling) or n=np (partly recaptured when downwelling)                    
            [n,np]=switchf(n,np);
            %commonterm=0.5*omega(np)*(W_ninp((2*NBIG-n)*LMU+1:(2*NBIG-n+1)*LMU,np)*WTMU)'; 
            commonterm=commonterm2;
            %% calculatin of Qm_kl,k--(i,n) and l-->(j,n')
            %% (1) downwelling to downwelling using PT               
            for m=0:M_Max               
                Pmat=[Fsr(np)*PT_MRAY(1:DStokes*LMU,1:DStokes*LMU,m+1)+Fsm(np)*PT_MMIE(1:DStokes*LMU,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+np-1))];
                Qmkl(DStokes*(2*NBIG-np)*LMU+1:DStokes*(2*NBIG-np+1)*LMU,DStokes*(2*NBIG-n)*LMU+1:DStokes*(2*NBIG-n+1)*LMU,m+1)=commonterm.*Pmat;                  


                %% (2) downwelling to upwelling using PR
                Pmat=[Fsr(np)*PR_MRAY(1:DStokes*LMU,1:DStokes*LMU,m+1)+Fsm(np)*PR_MMIE(1:DStokes*LMU,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+np-1))];
                Qmkl(DStokes*(np-1)*LMU+1:DStokes*np*LMU,DStokes*(2*NBIG-n)*LMU+1:DStokes*(2*NBIG-n+1)*LMU,m+1)=commonterm.*Pmat;                  
            end

            % switch back after job (2) is done
            [n,np]=switchf(n,np);  
            % job (2) done                        

        end % of np            
    end %of n


%Qmkl

    %% point 2 correction for adding method integration - for downwelling
    if NFIRST~=GROUPLAYER(NUMADD)  %correction if it is not the first adding             
        for n=1:NBIG
            SDF=WMU12;  
            for m=0:M_Max                  
                Qmkl(DStokes*(NBIG-1)*LMU+1:DStokes*NBIG*LMU,DStokes*(2*NBIG-n)*LMU+1:DStokes*(2*NBIG-n+1)*LMU,m+1)=SDF*W_ninp((2*NBIG-n)*DStokes*LMU+1:(2*NBIG-n+1)*DStokes*LMU,NBIG)'.*SLAST(1:DStokes*LMU,1:DStokes*LMU,m+1);      
            end
        end
    else
        for n=1:NBIG    

            SDF=2*XMU12.*WMU12;
            for m=0:M_Max  
                Qmkl(DStokes*(NBIG-1)*LMU+1:DStokes*NBIG*LMU,DStokes*(2*NBIG-n)*LMU+1:DStokes*(2*NBIG-n+1)*LMU,m+1)=SDF*W_ninp((2*NBIG-n)*DStokes*LMU+1:(2*NBIG-n+1)*DStokes*LMU,NBIG)'.*PRSur_M(1:DStokes*LMU,1:DStokes*LMU,m+1);                      
            end
        end            
    end    


    
    % Part 2: calculation of the R-matrix
    L=0;
    XME=XMU;
    XMJ=XMU;

    Rmle(1:DStokes*LXMUeALL,1:Qrown2,1:M_Max+1)=0;

    % Part A: upwelling    

    for np=1:NBIG          
        taunp=tau(np);
        taunp1=tau(np+1);
        dtnp=tau(np+1)-tau(np);

            for n=1:np % namely n<=np

                if n==np                    
                    FACT2A=exp(-taunp./XME12mat);
                    FACT2B=exp(-taunp1./XME12mat);

                    FACT1=0.25/dtnp./(XME12mat-XMJ12mat).*XMJ12mat;
                    FACT2=XME12mat./XMJ12mat.*(FACT2A-FACT2B);

                    FACT3A=exp(-tau(n)./XME12mat);
                    FACT3B=exp(-dtnp./XMJ12mat);
                    FACT3=FACT3A.*(1-FACT3B);  

                    SDF=omega(n)*FACT1.*(FACT2-FACT3);


                    % correction on uj=ue
                    for jj=1:LMU %incidence angle 
                        FACT1S=0.25/dtnp/XMJ(jj);
                        FACT2S=XMJ(jj)*exp(-taunp/XMJ(jj))*[1-exp(-dtnp/XMJ(jj))];
                        FACT3S=XMJ(jj)*exp(-taunp/XMJ(jj))*[dtnp/XMJ(jj)*exp(-dtnp/XMJ(jj))]; 

                        IDjj=DStokes*(jj-1)+1;
                        SDF(IDjj:IDjj+DStokes1,IDjj:IDjj+DStokes1)=omega(n)*FACT1S*(FACT2S-FACT3S);
                    end

                    for m=0:M_Max                          
                        pmat=[Fsr(n)*PT_MasRAY(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1)+Fsm(n)*PT_MasMIE(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+n-1))];
                        pSDF=SDF.*pmat;                                               
                        Rmle(1:DStokes*LXMUeALL,DStokes*(np-1)*LMU+1:DStokes*np*LMU,m+1)=Rmle(1:DStokes*LXMUeALL,DStokes*(np-1)*LMU+1:DStokes*np*LMU,m+1)+pSDF;
                    end

                else

                    NFACT1=0.25/dtnp./(XME12mat-XMJ12mat).*XMJ12mat;
                    NFACT2A=exp((tau(n+1)-taunp)./XMJ12mat-tau(n+1)./XME12mat);
                    NFACT2B=exp((tau(n+1)-taunp1)./XMJ12mat-tau(n+1)./XME12mat);
                    NFACT2=NFACT2A-NFACT2B;

                    NFACT3A=exp((tau(n)-taunp )./XMJ12mat-tau(n)./XME12mat);
                    NFACT3B=exp((tau(n)-taunp1)./XMJ12mat-tau(n)./XME12mat);
                    NFACT3=NFACT3A-NFACT3B; 


                    SDF=omega(n)*NFACT1.*(NFACT2-NFACT3);


                    for jj=1:LMU %incidence angle                            
                        NFACT1S=0.25/dtnp/XMJ(jj);
                        NFACT2S=tau(n+1)*(exp(-taunp/XMJ(jj))-exp(-taunp1/XMJ(jj)));
                        NFACT3S=tau(n)*(exp(-taunp/XMJ(jj))-exp(-taunp1/XMJ(jj))); 
                        IDjj=DStokes*(jj-1)+1; 
                        SDF(IDjj:IDjj+DStokes1,IDjj:IDjj+DStokes1)=omega(n)*NFACT1S*(NFACT2S-NFACT3S);
                    end

                    for m=0:M_Max                          
                        pmat=[Fsr(n)*PT_MasRAY(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1)+Fsm(n)*PT_MasMIE(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+n-1))];
                        pSDF=SDF.*pmat;
                        Rmle(1:DStokes*LXMUeALL,DStokes*(np-1)*LMU+1:DStokes*np*LMU,m+1)=Rmle(1:DStokes*LXMUeALL,DStokes*(np-1)*LMU+1:DStokes*np*LMU,m+1)+pSDF;
                    end

                end
            end        
    end %of np

%     % Part B: downwelling -- not necessary to calculate due to later surface-reflection effect correction



%% point 3 correction for adding method integration       
%if NFIRST~=GROUPLAYER(NUMADD)  %correction if it is not the first adding
    %% Part 1: Upwelling  
    Rmle(1:DStokes*LXMUeALL,DStokes*LMU*(NBIG-1)+1:DStokes*LMU*NBIG,1:M_Max+1)=0;


    for n=1:NBIG-1 %NFINAL
        FA=exp((tau(n+1)-tau(NBIG))./XMJ12mat-tau(n+1)./XME12mat);
        FB=exp((tau(n)  -tau(NBIG))./XMJ12mat-tau(n)  ./XME12mat);

        SDF=0.25*omega(n)./(XME12mat-XMJ12mat).*(FA-FB);
        % correction on uj=ue 
        for jj=1:LMU %incidence angle
            F0=0.25/XME(jj)/XMJ(jj)*omega(n)*exp(-tau(NBIG)/XMJ(jj));
            IDjj=DStokes*(jj-1)+1;
            SDF(IDjj:IDjj+DStokes1,IDjj:IDjj+DStokes1)=F0*(tau(n+1)-tau(n));
        end

        for m=0:M_Max                              
            Pmat=[Fsr(n)*PT_MasRAY(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1)+Fsm(n)*PT_MasMIE(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+n-1))];
            TRmle=SDF.*Pmat;
            Rmle(1:DStokes*LXMUeALL,DStokes*LMU*(NBIG-1)+1:DStokes*LMU*NBIG,m+1)=Rmle(1:DStokes*LXMUeALL,DStokes*LMU*(NBIG-1)+1:DStokes*LMU*NBIG,m+1)+TRmle;
        end

    end


%% Part 2: Downwelling (including the scattering properties from the last layer) 

    Rmle(1:DStokes*LXMUeALL,Qrown2+1:Qrown,1:M_Max+1)=0;
    for np=1:NBIG-1 %NFINAL
        taunp=tau(np);
        taunp1=tau(np+1);
        dtnp=tau(np+1)-tau(np);
         for n=np:NBIG-1 % namely n>=np
            if n==np                     
                FACT2A=exp(-taunp./XME12mat);
                FACT2B=exp(-taunp1./XME12mat);

                FACT1=0.25/dtnp./(XME12mat+XMJ12mat).*XMJ12mat;
                FACT2=XME12mat./XMJ12mat.*(FACT2A-FACT2B);

                FACT3A=exp(-tau(n+1)./XME12mat);
                FACT3B=exp(-dtnp./XMJ12mat);
                FACT3=FACT3A.*(1-FACT3B);  

                SDF=omega(n)*FACT1.*(FACT2-FACT3);
                for m=0:M_Max  
                    pmat=[Fsr(n)*PR_MRAY(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1)+Fsm(n)*PR_MMIE(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+n-1))];
                    pSDF=SDF.*pmat;
                    Rmle(1:DStokes*LXMUeALL,(2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,m+1)=Rmle(1:DStokes*LXMUeALL,(2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,m+1)+pSDF;
                end
            else                                   
                NFACT1=0.25/dtnp./(XME12mat+XMJ12mat).*XMJ12mat;

                NFACT2A=exp((taunp1-tau(n))./XMJ12mat-tau(n)./XME12mat);
                NFACT2B=exp((taunp-tau(n))./XMJ12mat-tau(n)./XME12mat);
                NFACT2=NFACT2A-NFACT2B;

                NFACT3A=exp((taunp1-tau(n+1))./XMJ12mat-tau(n+1)./XME12mat);
                NFACT3B=exp((taunp -tau(n+1))./XMJ12mat-tau(n+1)./XME12mat);
                NFACT3=NFACT3A-NFACT3B; 

                SDF=omega(n)*NFACT1.*(NFACT2-NFACT3); 
                for m=0:M_Max  
                    pmat=[Fsr(n)*PR_MRAY(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1)+Fsm(n)*PR_MMIE(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1,MiemixtypeforL(NFIRST+n-1))];
                    pSDF=SDF.*pmat;
                    Rmle(1:DStokes*LXMUeALL,(2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,m+1)=Rmle(1:DStokes*LXMUeALL,(2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,m+1)+pSDF;
                end
            end
         end  
            % contribution of the last layer (a reflecting surface)
            if NFIRST~=GROUPLAYER(NUMADD)
                Fmatcorr=0.5./XME12mat;

                F0=XMJ12mat/(tau(np+1)-tau(np)).*exp(-tau(NBIG)./XME12mat);
                FA=exp((tau(np+1)-tau(NBIG))./XMJ12mat);
                FB=exp((tau(np)-tau(NBIG))./XMJ12mat);
                SDF=F0.*(FA-FB);
                for m=0:M_Max                         
                    pmat=SLAST(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1);
                    pSDF=pmat.*SDF;
                    Rmle(1:DStokes*LXMUeALL,(2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,m+1)=Rmle(1:DStokes*LXMUeALL,(2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,m+1) + pSDF.*Fmatcorr;                        
                end
            else
                % addition of reflection from the last layer - surface

                F0=XMJ12mat/(tau(np+1)-tau(np)).*exp(-tau(NBIG)./XME12mat);
                FA=exp((tau(np+1)-tau(NBIG))./XMJ12mat);
                FB=exp((tau(np)-tau(NBIG))./XMJ12mat);
                SDF=F0.*(FA-FB);

                for m=0:M_Max                         
                    pmat=PRSur_M(1:DStokes*LXMUeALL,1:DStokes*LMU,m+1);
                    pSDF=pmat.*SDF;
                    Rmle(1:DStokes*LXMUeALL,(2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,m+1)=Rmle(1:DStokes*LXMUeALL,(2*NBIG-np)*DStokes*LMU+1:(2*NBIG-np+1)*DStokes*LMU,m+1) + pSDF;
                end

            end     

    end % of np
    %aa=onescosmatR.*dRda;           
    Rmle(1:DStokes*LXMUeALL,NBIG*DStokes*LMU+1:(NBIG+1)*DStokes*LMU,1:M_Max+1)=0.0;


    %end  %of "if NFIRST~=GROUPLAYER(NUMADD)"    


    % Part 3: calculation of the intial distribution
    for n=1:NBIG           
       % Part A: down to upwelling
        dtn=tau(n+1)-tau(n);              
        SDF=0.5*exp(-tau(n)./XMU0ALLmat).*XMU0ALLmat./(XMU0ALLmat+XMUIALLmat)*omega(n).*WMUIALLmat.*(1-exp(-dtn*(1./XMU0ALLmat + 1./XMUIALLmat)))./(1-exp(-dtn./XMUIALLmat))*dtn;
        for m=0:M_Max
            SOURCE1(DStokes*(n-1)*LMU+1:DStokes*n*LMU,1:DStokes*LXMU0ALL,m+1)=SDF.*[Fsr(n)*PR_MRAY(1:DStokes*LMU,:,m+1)+Fsm(n)*PR_MMIE(1:DStokes*LMU,:,m+1,MiemixtypeforL(NFIRST+n-1))];                        
        end             

        % Part B: down to downwelling               
        NDN=NBIG-n+1;
        dtn=tau(NDN+1)-tau(NDN);
        SDF=0.5*exp(-tau(NDN)./XMU0ALLmat).*XMU0ALLmat./(XMUIALLmat-XMU0ALLmat)*omega(NDN).*WMUIALLmat.*(exp(-dtn./XMUIALLmat)-exp(-dtn./XMU0ALLmat))./(1-exp(-dtn./XMUIALLmat))*dtn; 



        % correction on ui=u0
        for i0=1:LMU %incidence angle 
            IDi0=DStokes*(i0-1)+1;
            XMU0=XMJ(i0);
            SDF(IDi0:IDi0+DStokes1,IDi0:IDi0+DStokes1)=0.5*exp(-tau(NDN)/XMJ(i0))*omega(NDN)*WTMU(i0)*dtn/XMU0*exp(-dtn/XMU0)/(1-exp(-dtn/XMU0))*dtn;            

        end
       % ------------------- %    
        for m=0:M_Max                   
            SOURCE2(DStokes*(n-1)*LMU+1:DStokes*n*LMU,1:DStokes*LXMU0ALL,m+1)=SDF.*[Fsr(NDN)*PT_MRAY(1:DStokes*LMU,:,m+1)+Fsm(NDN)*PT_MMIE(1:DStokes*LMU,:,m+1,MiemixtypeforL(NFIRST+NDN-1))];                    
        end
    end
    SOURCE=[SOURCE1(1:Qrown/2,1:DStokes*LXMU0ALL,1:M_Max+1);SOURCE2(1:Qrown/2,1:DStokes*LXMU0ALL,1:M_Max+1)];      


        %% point 4 correction for adding method integration
    if NFIRST~=GROUPLAYER(NUMADD)  %correction if it is not the first adding
        for m=0:M_Max 
            SOURCE((NBIG-1)*DStokes*LMU+1:NBIG*DStokes*LMU,:,m+1)=exp(-tau(NBIG)./XMU0ALLmat).*XMU0ALLmat.*SLAST(1:DStokes*LMU,:,m+1).*WMUIALLmat;
        end                                                     
    else
        for m=0:M_Max                 
            SOURCE((NBIG-1)*DStokes*LMU+1:NBIG*DStokes*LMU,:,m+1)=exp(-tau(NBIG)./XMU0ALLmat).*XMU0ALLmat*2.*XMUIALLmat.*PRSur_M(1:DStokes*LMU,:,m+1).*WMUIALLmat;
        end

    end   

    for m=0:M_Max             
        if m==0                                         
            SOURCE(:,:,m+1)=SOURCE(1:Qrown,1:DStokes*LXMU0ALL,m+1);   
        else
            % seperation of cosine and sine mode                    
            SOURCEs(:,:,m+1)=onessinmatS(1:Qrown,1:DStokes*LXMU0ALL).*SOURCE(1:Qrown,1:DStokes*LXMU0ALL,m+1);
            SOURCE(:,:,m+1)=onescosmatS(1:Qrown,1:DStokes*LXMU0ALL).*SOURCE(1:Qrown,1:DStokes*LXMU0ALL,m+1);   
        end
    end




    % seperation of cosine and sine mode
    for m=0:M_Max
        if m==0
           Rle(1:DStokes*LXMUeALL,1:Qrown,m+1)=Rmle(1:DStokes*LXMUeALL,1:Qrown,m+1);
           ADD=SOURCE(:,:,m+1);
           QPI(1:Qrown,1:DStokes*LXMU0ALL,m+1)=ADD;
           for np=1:NMultiplication
               ADD=Qmkl(1:Qrown,1:Qrown,m+1)*ADD;
               QPI(1:Qrown,1:DStokes*LXMU0ALL,m+1)=QPI(1:Qrown,1:DStokes*LXMU0ALL,m+1)+ADD;
           end           
        else
           Rles(1:DStokes*LXMUeALL,1:Qrown,m+1)=onessinmatR(1:DStokes*LXMUeALL,1:Qrown).*Rmle(1:DStokes*LXMUeALL,1:Qrown,m+1);
           Rle(1:DStokes*LXMUeALL,1:Qrown,m+1)=onescosmatR(1:DStokes*LXMUeALL,1:Qrown).*Rmle(1:DStokes*LXMUeALL,1:Qrown,m+1);         
           Qkls=onessinmatQ(1:Qrown,1:Qrown).*Qmkl(1:Qrown,1:Qrown,m+1);
           ADD=SOURCE(:,:,m+1)+SOURCEs(:,:,m+1);
           QPI(:,:,m+1)=ADD;
           for np=1:NMultiplication
               ADD=Qmkl(1:Qrown,1:Qrown,m+1)*ADD(1:Qrown,1:DStokes*LXMU0ALL)-2*Qkls*[onessinmatS.*ADD(1:Qrown,1:DStokes*LXMU0ALL)];
               QPI(:,:,m+1)=QPI(:,:,m+1)+ADD;               
           end           
        end

    end



        %% point 5 correction for adding method integration     
    if NFIRST~=GROUPLAYER(NUMADD)
        commat=-(1./XMU0OUT+1./XMUiOUT);
        EF=exp(tau(NBIG)*commat);
        for m=0:M_Max 
            if m==0
                tc=Rle(1:DStokes*LXMUeALL,1:Qrown,m+1)*QPI(1:Qrown,1:DStokes*LXMU0ALL,m+1); 
                Imc(:,:,m+1)=Imc(:,:,m+1).*EF+tc;                
                for jjj=1:Lfaipfai0
                    Im(:,1:LXMU0ALL,jjj,m+1)=Im(:,1:LXMU0ALL,jjj,m+1).*EF(:,1:DStokes:DStokes*LXMU0ALL)+tc(:,1:DStokes:DStokes*LXMU0ALL)*cos(m*faipfai0(jjj)); 
                end

            else           
                tcs=Rmle(:,:,m+1)*QPI(:,:,m+1)-2*Rles(:,:,m+1)*[onessinmatS.*QPI(1:Qrown,1:DStokes*LXMU0ALL,m+1)];               
                tc(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL)=onescosmatRS.*tcs;
                ts(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL)=onessinmatRS.*tcs;                 
                
                % !!! Following operation has to be put after the derivative calculation
                Imc(:,:,m+1)=Imc(:,:,m+1).*EF+tc; 
                Ims(:,:,m+1)=Ims(:,:,m+1).*EF+ts;                  

                for jjj=1:Lfaipfai0                
                    Im(:,1:LXMU0ALL,jjj,m+1)=Im(:,1:LXMU0ALL,jjj,m+1).*EF(:,1:DStokes:DStokes*LXMU0ALL)+[tc(:,1:DStokes:DStokes*LXMU0ALL)*cos(m*faipfai0(jjj))+ts(:,1:DStokes:DStokes*LXMU0ALL)*sin(m*faipfai0(jjj))]; 
                end
            end 
        end
    else   
        for m=0:M_Max                
            if m==0
                %Imc(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL,m+1)=Xm(1:DStokes*LXMUeALL,1:Qrown,m+1)*SOURCE(:,:,m+1);
                Imc(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL,m+1)=Rle(1:DStokes*LXMUeALL,1:Qrown,m+1)*QPI(1:Qrown,1:DStokes*LXMU0ALL,m+1);
                for jjj=1:Lfaipfai0
                    Im(:,1:LXMU0ALL,jjj,m+1)=Imc(:,1:DStokes:DStokes*LXMU0ALL,m+1)*cos(m*faipfai0(jjj));
                end
            else
                %Imcs=[Rle(:,:,m+1),-Rles(:,:,m+1);Rles(:,:,m+1),Rle(:,:,m+1)]*QPI(:,:,m+1);
                Imcs=Rmle(:,:,m+1)*QPI(:,:,m+1)-2*Rles(:,:,m+1)*[onessinmatS.*QPI(1:Qrown,1:DStokes*LXMU0ALL,m+1)];               
                Imc(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL,m+1)=onescosmatRS.*Imcs;
                Ims(1:DStokes*LXMUeALL,1:DStokes*LXMU0ALL,m+1)=onessinmatRS.*Imcs; 
              
                for jjj=1:Lfaipfai0
                    Im(:,1:LXMU0ALL,jjj,m+1)=Imc(:,1:DStokes:DStokes*LXMU0ALL,m+1)*cos(m*faipfai0(jjj))+Ims(:,1:DStokes:DStokes*LXMU0ALL,m+1)*sin(m*faipfai0(jjj));                     
                end
            end
        end
    end

    matcom=2*XMUeOUT./XMU0OUT; %change 39
    for m=0:M_Max         
        if m==0          
            SLAST(:,:,m+1)=Imc(:,:,m+1).*matcom;

        else
            SLAST(:,:,m+1)=[Imc(:,:,m+1)+Ims(:,:,m+1)].*matcom;
        end
    end
    %plot(a(:,i0))
    %hold on



    %% point 6 correction for adding method integration
    %%  ---  THIS BLOCK PRESERVES THE RESULTS OF THE PRESENT CALCULATIO FOR THE NEXT "ADDING" 
    %%  ---  THE TRANSITION PROBABILITIES INCLUDE BOTH SINGLE AND HIGHER ORDER SCATTERING.
    if NFIRST~=1 % Xu added: when NFIRST the whole layer is no longer necessary to compress into one single layer since final result has been obtained.
                 % This is where Larry's code crashes.
        % addition of single scattering in atmospheric layers
        for n=NFIRST:GROUPLAYER(end)-1 %NUMBIG
            taun=-[TAUALL(NFIRST-1)-TAUALL(n-1)];
            taun1=-[TAUALL(NFIRST-1)-TAUALL(n)];        


            ADF=[exp(-taun*(1./XMUeOUT+1./XMU0OUT))-exp(-taun1*(1./XMUeOUT+1./XMU0OUT))];
            SLASTdp=0.5*OMEGA(n)*XMUeOUT./(XMU0OUT+XMUeOUT).*ADF;

            for m=0:M_Max                 
                pmat=[RAYF(n)*PR_MRAY(:,:,m+1)+MIEF(n)*PR_MMIE(:,:,m+1,MiemixtypeforL(n))];
                SLAST(:,:,m+1)=SLAST(:,:,m+1)+SLASTdp.*pmat;
            end                          

        end               
    else
        % addition of single scattering in atmospheric layers
        for n=NFIRST:GROUPLAYER(end)-1 %NUMBIG
            if n~=1
                taun=-[0-TAUALL(n-1)];
                taun1=-[0-TAUALL(n)];
            else
                taun=-[0-0];
                taun1=-[0-TAUALL(n)];
            end                        
            ADF=[exp(-taun*(1./XMUeOUT+1./XMU0OUT))-exp(-taun1*(1./XMUeOUT+1./XMU0OUT))];
            SLASTdp=0.5*OMEGA(n)*XMUeOUT./(XMU0OUT+XMUeOUT).*ADF;


            for m=0:M_Max 
                pmat=[RAYF(n)*PR_MRAY(:,:,m+1)+MIEF(n)*PR_MMIE(:,:,m+1,MiemixtypeforL(n))];
                SLAST(:,:,m+1)=SLAST(:,:,m+1)+SLASTdp.*pmat;
            end

        end             
    end     

    % addition of single scattering via relfection from surface
    tau0=taun1;       
    SDF=2*XMUeOUT.*exp(-tau0*(1./XMUeOUT+1./XMU0OUT));    
    for m=0:M_Max 
        SLAST(:,:,m+1)=SLAST(:,:,m+1)+SDF.*PRSur_M(:,:,m+1);           
    end            

    
    NLAST=NFIRST;    
    Rle=[];
    Rles=[];    
    SOURCE=[];
    SOURCEs=[];
    QPI=[];
    Rmle=[];

end
toc
        
for m=0:M_Max         
    if m==0
       deltam=1;
    else
       deltam=2;
    end
    for jjj=1:Lfaipfai0
%         for i0=1:LXMU0ALL 
%             IM(:,i0,jjj)=IM(:,i0,jjj)+deltam*Im(:,i0,jjj,m+1);                     
%         end
        IM(:,1:LXMU0ALL,jjj)=IM(:,1:LXMU0ALL,jjj)+deltam*Im(:,1:LXMU0ALL,jjj,m+1); 
    end        
end




tic
  

% Following is a copy from HomogeneousAtmosphereV_BRDF3
% single scattering calculation directly using original phase function
% phase function and its expansion series needs to be further changed for vertically inhomogeneous atmosphere
IS(1:DStokes*LXMUeALL,1:LXMU0ALL,1:Lfaipfai0)=0;

if TAUALL(1)~=0
    TAUALL=[0,TAUALL];
end
tau0=TAUALL(end-1);
NBIGALL=length(TAUALL)-2;


for jjj=1:Lfaipfai0
    for i0=1:LXMU0ALL %initial incidence angle     
        XMUI=XMU0ALL(i0);XMUJ=XMUeALL;
        XMVI=sqrt(1-XMUI^2);XMVJ=sqrt(1-XMUJ.^2);
        XIXJ=XMUI*XMUJ;VIVJ=XMVI*XMVJ;
        CSTHR=-XIXJ+VIVJ*cos(faipfai0(jjj));
        
        
        CSTHR(find(CSTHR==1))=0.999999999;
        CSTHR(find(CSTHR==-1))=-0.999999999;    
   


        % note that in following "Pij_nr", Rayleigh scattering contribution has been integrated.
        % Rayleigh part
        PR11=interp1(cos(theta),P11_nr,CSTHR,'spline');       
        PR21=interp1(cos(theta),P21_nr,CSTHR,'spline');  
        PR33=interp1(cos(theta),P33_nr,CSTHR,'spline');  
        PR43=interp1(cos(theta),P43_nr,CSTHR,'spline');          
        PR44=interp1(cos(theta),P44_nr,CSTHR,'spline');  %added 

        
% cosi1R=[XMUJ*XMVI+XMUI*XMVJ.*cos(faipfai0(jjj))]./sqrt(1-CSTHR.^2);        
% cosi2R=[-XMUI*XMVJ-XMUJ*XMVI.*cos(faipfai0(jjj))]./sqrt(1-CSTHR.^2); 

        cosi1R=[XMUJ+XMUI*CSTHR]./sqrt(1-CSTHR.^2)/XMVI;        
        cosi2R=[-XMUI-XMUJ.*CSTHR]./sqrt(1-CSTHR.^2)./XMVJ; 
        sini1R=[XMVJ.*sin(-faipfai0(jjj))]./sqrt(1-CSTHR.^2); 
        sini2R=[XMVI*sin(-faipfai0(jjj))]./sqrt(1-CSTHR.^2);             
        cos2alfa0=2*cosi1R.^2-1;
        sin2alfa0=2*sini1R.*cosi1R;

        cos2alfa=2*cosi2R.^2-1;
        sin2alfa=2*sini2R.*cosi2R;          
        
        [P11,P12,P13,P21,P22,P23,P24,P31,P32,P33,P34,P42,P43,P44]=RPR_Cal(PR11,PR21,PR11,PR33,-PR43,PR44,cos2alfa0,sin2alfa0,cos2alfa,sin2alfa);        
        
        % Mie part
        
        for L=1:LMiemixtype
            PR11m=interp1(cos(theta),P11_nm(L,:),CSTHR,'spline');       
            PR21m=interp1(cos(theta),P21_nm(L,:),CSTHR,'spline');  
            PR33m=interp1(cos(theta),P33_nm(L,:),CSTHR,'spline');  
            PR43m=interp1(cos(theta),P43_nm(L,:),CSTHR,'spline');          
            PR44m=interp1(cos(theta),P44_nm(L,:),CSTHR,'spline');  %added        
            [P11m(L,:),P12m(L,:),P13m(L,:),...
             P21m(L,:),P22m(L,:),P23m(L,:),P24m(L,:),...
             P31m(L,:),P32m(L,:),P33m(L,:),P34m(L,:),...
             P42m(L,:),P43m(L,:),P44m(L,:)]=RPR_Cal(PR11m,PR21m,PR11m,PR33m,-PR43m,PR44m,cos2alfa0,sin2alfa0,cos2alfa,sin2alfa);        

        end        
        
        
        [P11sur,P12sur,P22sur,P33sur,P34sur,P44sur]=Pmat_PolaBRDF_SurfaceRPV_Liz_OptOa(nn,vv,epsirol,r01,b1,k1,XMUI,XMUJ,cos(faipfai0(jjj)),LMU);      
      
        [RPRsur11,RPRsur12,RPRsur13,...
         RPRsur21,RPRsur22,RPRsur23,RPRsur24,...
         RPRsur31,RPRsur32,RPRsur33,RPRsur34,...
                  RPRsur42,RPRsur43,RPRsur44]=RPR_Cal(P11sur,P12sur,P22sur,P33sur,P34sur,P44sur,cos2alfa0,sin2alfa0,cos2alfa,sin2alfa);        

       

              
              
%         if Lambertian==1 % lambertian surface
%             CF=pi;
%         else
%             CF=pi./XMU;
%         end
        CF=pi./XMUeALL;
        UUI=1./XMUeALL+1/XMUI;
        EC=exp(-tau0*UUI);
        XEF=XMUI*EC.*CF;
        PRSur(1:DStokes:DStokes*LXMUeALL,i0,jjj)=XEF.*RPRsur11;
        PRSur(2:DStokes:DStokes*LXMUeALL,i0,jjj)=XEF.*RPRsur21;
        PRSur(3:DStokes:DStokes*LXMUeALL,i0,jjj)=XEF.*RPRsur31;

        
        for n=1:NBIGALL           
            
            IStmp0=OMEGA(n)/4*XMUI./(XMUeALL+XMUI).*[exp(-TAUALL(n)*UUI)-exp(-TAUALL(n+1)*UUI)]; 
            IStmp11=IStmp0.*[RAYF(n)*P11+MIEF(n)*P11m(MiemixtypeforL(n),:)];           
            IStmp21=IStmp0.*[RAYF(n)*P21+MIEF(n)*P21m(MiemixtypeforL(n),:)];  
            IStmp31=IStmp0.*[RAYF(n)*P31+MIEF(n)*P31m(MiemixtypeforL(n),:)];
            
            
            IS(1:DStokes:DStokes*LXMUeALL,i0,jjj)=IS(1:DStokes:DStokes*LXMUeALL,i0,jjj)+IStmp11'; 
            IS(2:DStokes:DStokes*LXMUeALL,i0,jjj)=IS(2:DStokes:DStokes*LXMUeALL,i0,jjj)+IStmp21';
            IS(3:DStokes:DStokes*LXMUeALL,i0,jjj)=IS(3:DStokes:DStokes*LXMUeALL,i0,jjj)+IStmp31';
       
        end
        
        
        if DStokes==4
            PRSur(4:DStokes:DStokes*LXMUeALL,i0,jjj)=0;            
            for n=1:NBIGALL           
                IStmp41=0;
                IS(4:DStokes:DStokes*LXMUeALL,i0,jjj)=IS(4:DStokes:DStokes*LXMUeALL,i0,jjj)+IStmp41';
            end                        
        end
        
        

        

    end
end % of jjj

IS=IS+PRSur;%
toc
