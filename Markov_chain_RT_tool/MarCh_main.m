%% The program is developed to sumilate top-of-atmosphere polarized radiance  
%% in a coupled atmosphere-surface system, with application to  
%% 1) Multi-angle Imaging SpectroRadiometer
%% 2) Airborne Multi-angle SpectroPolarimetric Imager (AirMSPI)
%% 3) Multi-Angle Imager for Aerosols

%% Last modified 07/17/2017

%% The program is based on the following references
%% 1. Xu et al., Markov chain formalism for vector radiative transfer in a plane-parallel atmosphere overlying a polarizing surface, Opt. Lett. 36, 2083-2085, 2011.
%% 2. Xu et al., Linearization of Markov chain formalism for vector radiative transfer in a plane-parallel atmosphere/surface system, Appl. Opt. 51, 3491-3507, 2012.
%% 3. Xu et al., Markov chain formalism for polarized light transfer in plane-parallel atmospheres, with numerical comparison to the Monte Carlo method, Opt. Express 19, 946-967, 2011.

%% The code was developed by Feng Xu (F.Xu, FENG.XU@JPL.NASA.GOV) at the Jet Propulsion Laboratory, Pasadena.  
 
%% Notes: 
%% 1) Bugs, revisions and/or improvements of the code should be reported to the author (F. Xu)
%% 2) In any publication or operational retrieval algorithm using the code, the source of the code is requested to be acknowledged with proper references being made.
%% 3) Code distribution beyond MISR/AirMSPI/MAIA group use needs approval from the author

clear
clear global SSA_p_TYPE Aerosol_Type MiemixtypeforL Miemixtype LMiemixtype AF_p theta TAUALL RAYF OMEGA P11_nm_TYPE P21_nm_TYPE P22_nm_TYPE P33_nm_TYPE P43_nm_TYPE P44_nm_TYPE
global Aerosol_Type SSA_p_TYPE MiemixtypeforL Miemixtype LMiemixtype AF_p theta TAUALL RAYF OMEGA P11_nm_TYPE P21_nm_TYPE P22_nm_TYPE P33_nm_TYPE P43_nm_TYPE P44_nm_TYPE
addpath(genpath('/Users/huikyole/Work/MAIA_UQ/YAMLMatlab_0.4.3'))
%% Input file for HomogeneousAtmosphereV_BRDF_Bob
choice=2;%1 for HG phase function described by (g1, g2) which can be different for different layer
%choice=3;%1 for HG phase function described by (g1, g2) which can be different for different layer
         %2 for Mie phase function described by (mr, mi, rm, lnS, tau_AOD) which can be different for different layer
         %3 for direct input of phase matrix
delta_use=1;% Turn on (1) or off (0) the delta-truncation for aerosol phase matrix

%% Part 2: Rayleigh and aerosol optical depth
%%      Layer Number           Rayleigh Optical Depth     Aerosol (TYPE 1) OD      Aerosol (TYPE 2) OD
OD = csvread('Inputs/OD.csv')
OD(:,4)=OD(:,3)*0.2;
OD(:,3)=OD(:,3)*0.8;
Aerosol_Type=size(OD,2)-2;

%% Part 2: Single scattering albedo and phase matrix
Ltheta=2001;
theta=linspace(0,pi,Ltheta); % theta for interpolation
costheta=cos(theta);
if choice==1 %% HG phase function described by (g1, g2) which can be different for different layer
    
    SSA_p_TYPE=[0.838,0.838]; % single scattering albedo 
    % Asymmetry factor of HG function: -1<g<1 (g > 0, forward scattering is dominant, g < 0, backscattering predominates.)
    g(1,:)=[0.5,0.7]; % Aerosol of type 1
    g(2,:)=[0.5,0.7]; % Aerosol of type 2
    
    f1=[0.143,0.143];
    h=[1,1]; %% Parameters defined by Eqs.(46)-(49) of C. J. Braak et al, JQSRT 69, 585(2001).
    p=[1,1];

    
    for n=1:Aerosol_Type
        P11R_a=HG_Phase(g(n,1),costheta);
        P11R_b=HG_Phase(g(n,2),costheta);
        P11_nm_TYPE(n,:)=f1(n)*P11R_a+(1-f1(n))*P11R_b;    

        % The phase matrix for intestellar particles are constructed according to Eqs.(46)-(49) of C. J. Braak et al, JQSRT 69, 585(2001).
        a2overa1(n,:)=h(n)+(1-h(n))*costheta;
        a3overa1(n,:)=2*costheta./(1+costheta.^2); %Note here correction on Eq.(46) was made according to the statement below Eq.(49) 
        b1overa1(n,:)=-p(n)*(1-costheta.^2)./(1+costheta.^2);
        b2overa1(n,1:Ltheta)=0;
        a4overa1(n,1:Ltheta)=0;  

        P22_nm_TYPE(n,:)=a2overa1(n,:).*P11_nm_TYPE(n,:);
        P33_nm_TYPE(n,:)=a3overa1(n,:).*P11_nm_TYPE(n,:);
        P21_nm_TYPE(n,:)=b1overa1(n,:).*P11_nm_TYPE(n,:);
        P43_nm_TYPE(n,:)=-b2overa1(n,1:Ltheta).*P11_nm_TYPE(n,:);
        P44_nm_TYPE(n,:)=a4overa1(n,1:Ltheta).*P11_nm_TYPE(n,:);  

    end
     
    
end

if choice==2 %% Mie phase function described by (mr, mi, rm, lnS, tau_AOD) which can be different for different layer
    %% refer to code "LM_VRT_Surface_2mode" & "V_MISR_Kokhanovsky_IQU2mode_Imp2b"    
    parameters=ReadYaml('choice2.yaml')
    lamda0=parameters.lamda0
    %lamda0=0.443e-6;     % incident wavelength
    m_r=[1.3800,1.4500]; % Real part of aerosol refractive index
    m_i=[0.1,0];         % Imaginary part of aerosol refractive index
    rm=[0.10e-6,0.10e-6];% Medium radius of Log-normal aerosol size distribution (ASD) 
    lnS=[1.0,1.0];       % Wiss of Log-normal ASD
    
    % First integral range over particle size in [r_min, r_max]
    r_min=0.05e-6;
    r_max=20.0e-6; % Maxmum aerosl size for cutting the ASD tail
    r_number=20;   % Number of intervals for Gauss integral in [r_min, r_max]
    gamma_dis=0;   % Set "gamma_dis=0" for Log-normal distribution

    % Gauss points for numerical integral
    ddr=(r_max-r_min)/r_number; % r_number is interval number
    g_number=100;% Gaussian point in each intercal
    LR=r_number*g_number;
    r(LR)=0;wtr(LR)=0;    

    for n=1:r_number     
        n_start=(n-1)*g_number+1;
        n_end=n*g_number;

        r_minp=r_min+ddr*(n-1);
        r_maxp=r_min+ddr*n;        
        [r0,wtr0]=lgwt(g_number,r_minp,r_maxp);
        r(n_start:n_end)=r0;
        wtr(n_start:n_end)=wtr0;
    end   
    % Second integral range over particle size in [0, r_min]    
    g_number0=30; % Number of Gauss points in [0, r_min]
    [r00,wtr00]=lgwt(g_number0,0,r_min);
    r=[r00,r];
    wtr=[wtr00,wtr];
    LR=LR+g_number0;


    x=(pi*2*r_max/lamda0);
    Na=1;
    Nb=4.05;
    Nc=0.34;
    Nd=8;
    N1=ceil(Na*x+Nb*x^Nc+Nd);

    paiangle(N1,Ltheta)=0;
    tauangle(N1,Ltheta)=0;    

    pai0=0;
    n=1;
    paiangle(1,1:Ltheta)=1;
    pai1=paiangle(1,1:Ltheta);   
    tauangle(1,1:Ltheta)=n*cos(theta).*pai1-(n+1)*pai0;% Eq.(9) of ref.[4]

    % derivatives
    n=2;
    paiangle(2,1:Ltheta)=(2*n-1)/(n-1)*cos(theta).*pai1-n/(n-1)*pai0;
    pai2=paiangle(2,1:Ltheta);      
    tauangle(2,1:Ltheta)=n*cos(theta).*pai2-(n+1)*pai1;        
    for n=3:N1
        paiangle(n,1:Ltheta)=(2*n-1)/(n-1)*cos(theta).*pai2-n/(n-1)*pai1;
        tauangle(n,1:Ltheta)=n*cos(theta).*paiangle(n,1:Ltheta)-(n+1)*pai2;  

        pai1=pai2;
        pai2=paiangle(n,1:Ltheta);    
    end 
    
    P11_nm(Aerosol_Type,Ltheta)=0;
    P21_nm(Aerosol_Type,Ltheta)=0;
    P22_nm(Aerosol_Type,Ltheta)=0;    
    P33_nm(Aerosol_Type,Ltheta)=0;
    P43_nm(Aerosol_Type,Ltheta)=0;
    P44_nm(Aerosol_Type,Ltheta)=0;

    
    for n=1:Aerosol_Type
        
        [P11_nm_TYPE(n,:),P21_nm_TYPE(n,:),P33_nm_TYPE(n,:),P43_nm_TYPE(n,:),P44_nm_TYPE(n,:),SSA_p_TYPE(n)]=Phase_mat_cal(lamda0,m_r(n),m_i(n),rm(n),lnS(n),r,wtr,LR,gamma_dis,paiangle,tauangle,Ltheta,N1,Na,Nb,Nc,Nd);  
        
        P22_nm_TYPE(n,:)=P11_nm_TYPE(n,:); % for Mie spherical particle                 
%         sumP11=0;
%         for thetan=1:length(theta)-1
%             sumP11=sumP11+P11_nm_TYPE(n,thetan)*[theta(thetan+1)-theta(thetan)]*sin(theta(thetan));
%         end
%         % check whether the integral result is euqal to 2 (apprixmately)
%         sumP11         
    end
    
end

if choice==3
    SSA_p_TYPE=[0.838,0.838]; % single scattering albedo   
    %% Eaxmple: Phase Matrix from Olga for nonspherical particles (Provided by Olga)
    %  ANGLE         F11         F21          F22          F33          F43          F44
    Pmat(:,:,1)=[      0   7.5230e+02            0   7.5226e+02   7.5226e+02            0   7.5223e+02
    1.7453e-02   5.3903e+02   4.5521e-02   5.3900e+02   5.3899e+02   2.0724e-01   5.3896e+02
    3.4907e-02   3.1063e+02   7.5337e-02   3.1060e+02   3.1059e+02   6.3844e-01   3.1056e+02
    5.2360e-02   1.7953e+02   7.0600e-02   1.7950e+02   1.7949e+02   9.6488e-01   1.7946e+02
    6.9813e-02   1.0769e+02   6.3074e-02   1.0766e+02   1.0764e+02   1.0634e+00   1.0762e+02
    8.7266e-02   6.7586e+01   6.0907e-02   6.7562e+01   6.7533e+01   1.0106e+00   6.7514e+01
    1.0472e-01   4.4388e+01   6.2463e-02   4.4365e+01   4.4330e+01   8.9877e-01   4.4314e+01
    1.2217e-01   3.0441e+01   6.6404e-02   3.0420e+01   3.0380e+01   7.7113e-01   3.0367e+01
    1.3963e-01   2.1770e+01   7.1205e-02   2.1749e+01   2.1707e+01   6.4744e-01   2.1696e+01
    1.5708e-01   1.6198e+01   7.5782e-02   1.6177e+01   1.6134e+01   5.3651e-01   1.6125e+01
    1.7453e-01   1.2500e+01   7.9551e-02   1.2479e+01   1.2435e+01   4.4029e-01   1.2428e+01
    1.9199e-01   9.9735e+00   8.2159e-02   9.9521e+00   9.9086e+00   3.5876e-01   9.9023e+00
    2.0944e-01   8.1998e+00   8.3704e-02   8.1777e+00   8.1349e+00   2.9108e-01   8.1298e+00
    2.2689e-01   6.9828e+00   8.3787e-02   6.9594e+00   6.9166e+00   2.4144e-01   6.9126e+00
    2.4435e-01   5.9780e+00   8.3877e-02   5.9542e+00   5.9133e+00   1.9178e-01   5.9103e+00
    2.6180e-01   5.2587e+00   8.2667e-02   5.2339e+00   5.1940e+00   1.5609e-01   5.1919e+00
    2.7925e-01   4.6960e+00   8.0724e-02   4.6701e+00   4.6311e+00   1.2685e-01   4.6299e+00
    2.9671e-01   4.2433e+00   7.8056e-02   4.2163e+00   4.1783e+00   1.0312e-01   4.1778e+00
    3.1416e-01   3.8677e+00   7.4720e-02   3.8396e+00   3.8023e+00   8.3871e-02   3.8027e+00
    3.3161e-01   3.5485e+00   7.0892e-02   3.5192e+00   3.4826e+00   6.7936e-02   3.4838e+00
    3.4907e-01   3.2729e+00   6.6823e-02   3.2424e+00   3.2064e+00   5.4484e-02   3.2084e+00
    3.6652e-01   3.0326e+00   6.2796e-02   3.0010e+00   2.9654e+00   4.3039e-02   2.9681e+00
    3.8397e-01   2.8220e+00   5.9073e-02   2.7892e+00   2.7539e+00   3.3353e-02   2.7574e+00
    4.0143e-01   2.6360e+00   5.5725e-02   2.6020e+00   2.5669e+00   2.5225e-02   2.5711e+00
    4.1888e-01   2.4699e+00   5.2703e-02   2.4348e+00   2.3997e+00   1.8389e-02   2.4046e+00
    4.3633e-01   2.3198e+00   4.9880e-02   2.2836e+00   2.2484e+00   1.2570e-02   2.2540e+00
    4.5379e-01   2.1827e+00   4.7177e-02   2.1454e+00   2.1101e+00   7.5430e-03   2.1164e+00
    4.7124e-01   2.0562e+00   4.4572e-02   2.0179e+00   1.9823e+00   3.1503e-03   1.9893e+00
    4.8869e-01   1.9388e+00   4.2093e-02   1.8995e+00   1.8636e+00  -7.0598e-04   1.8712e+00
    5.0615e-01   1.8295e+00   3.9784e-02   1.7892e+00   1.7530e+00  -4.1396e-03   1.7612e+00
    5.2360e-01   1.7274e+00   3.7680e-02   1.6862e+00   1.6495e+00  -7.2578e-03   1.6584e+00
    5.4105e-01   1.6320e+00   3.5778e-02   1.5899e+00   1.5528e+00  -1.0072e-02   1.5622e+00
    5.5851e-01   1.5427e+00   3.4043e-02   1.4998e+00   1.4622e+00  -1.2577e-02   1.4722e+00
    5.7596e-01   1.4588e+00   3.2417e-02   1.4152e+00   1.3770e+00  -1.4804e-02   1.3876e+00
    5.9341e-01   1.3799e+00   3.0906e-02   1.3356e+00   1.2968e+00  -1.6775e-02   1.3079e+00
    6.1087e-01   1.3054e+00   2.9462e-02   1.2604e+00   1.2211e+00  -1.8511e-02   1.2327e+00
    6.2832e-01   1.2350e+00   2.8070e-02   1.1895e+00   1.1495e+00  -2.0056e-02   1.1615e+00
    6.4577e-01   1.1684e+00   2.6739e-02   1.1223e+00   1.0817e+00  -2.1478e-02   1.0942e+00
    6.6323e-01   1.1056e+00   2.5526e-02   1.0591e+00   1.0177e+00  -2.2849e-02   1.0307e+00
    6.8068e-01   1.0465e+00   2.4451e-02   9.9953e-01   9.5527e-01  -2.4225e-02   9.7081e-01
    6.9813e-01   9.9096e-01   2.3573e-02   9.4362e-01   8.9701e-01  -2.5593e-02   9.1454e-01
    7.1558e-01   9.3892e-01   2.2898e-02   8.9127e-01   8.4280e-01  -2.6898e-02   8.6178e-01
    7.3304e-01   8.9013e-01   2.2372e-02   8.4221e-01   7.9233e-01  -2.8109e-02   8.1229e-01
    7.5049e-01   8.4433e-01   2.1964e-02   7.9620e-01   7.4528e-01  -2.9227e-02   7.6580e-01
    7.6794e-01   8.0129e-01   2.1660e-02   7.5299e-01   7.0138e-01  -3.0267e-02   7.2207e-01
    7.8540e-01   7.6079e-01   2.1431e-02   7.1236e-01   6.6037e-01  -3.1242e-02   6.8091e-01
    8.0285e-01   7.2266e-01   2.1272e-02   6.7414e-01   6.2203e-01  -3.2173e-02   6.4212e-01
    8.2030e-01   6.8681e-01   2.1194e-02   6.3824e-01   5.8623e-01  -3.3073e-02   6.0561e-01
    8.3776e-01   6.5310e-01   2.1173e-02   6.0451e-01   5.5280e-01  -3.3943e-02   5.7125e-01
    8.5521e-01   6.2142e-01   2.1200e-02   5.7284e-01   5.2159e-01  -3.4772e-02   5.3758e-01
    8.7266e-01   5.9165e-01   2.1273e-02   5.4311e-01   4.9245e-01  -3.5550e-02   5.0621e-01
    8.9012e-01   5.6369e-01   2.1391e-02   5.1521e-01   4.6375e-01  -3.6283e-02   4.7700e-01
    9.0757e-01   5.3742e-01   2.1527e-02   4.8903e-01   4.3677e-01  -3.6994e-02   4.4978e-01
    9.2502e-01   5.1276e-01   2.1690e-02   4.6447e-01   4.1141e-01  -3.7704e-02   4.2443e-01
    9.4248e-01   4.8963e-01   2.1882e-02   4.4147e-01   3.8760e-01  -3.8420e-02   4.0084e-01
    9.5993e-01   4.6796e-01   2.2115e-02   4.1994e-01   3.6527e-01  -3.9136e-02   3.7890e-01
    9.7738e-01   4.4768e-01   2.2384e-02   3.9981e-01   3.4434e-01  -3.9846e-02   3.5850e-01
    9.9484e-01   4.2870e-01   2.2670e-02   3.8099e-01   3.2473e-01  -4.0540e-02   3.3953e-01
    1.0123e+00   4.1092e-01   2.2972e-02   3.6339e-01   3.0633e-01  -4.1219e-02   3.2188e-01
    1.0297e+00   3.9427e-01   2.3279e-02   3.4692e-01   2.8906e-01  -4.1895e-02   3.0545e-01
    1.0472e+00   3.7869e-01   2.3597e-02   3.3152e-01   2.7288e-01  -4.2565e-02   2.9016e-01
    1.0647e+00   3.6411e-01   2.3933e-02   3.1713e-01   2.5769e-01  -4.3227e-02   2.7508e-01
    1.0821e+00   3.5045e-01   2.4275e-02   3.0365e-01   2.4344e-01  -4.3876e-02   2.6094e-01
    1.0996e+00   3.3767e-01   2.4623e-02   2.9105e-01   2.3008e-01  -4.4508e-02   2.4769e-01
    1.1170e+00   3.2571e-01   2.4977e-02   2.7928e-01   2.1756e-01  -4.5117e-02   2.3528e-01
    1.1345e+00   3.1450e-01   2.5329e-02   2.6826e-01   2.0582e-01  -4.5697e-02   2.2364e-01
    1.1519e+00   3.0400e-01   2.5697e-02   2.5795e-01   1.9478e-01  -4.6238e-02   2.1272e-01
    1.1694e+00   2.9416e-01   2.6061e-02   2.4830e-01   1.8442e-01  -4.6739e-02   2.0249e-01
    1.1868e+00   2.8492e-01   2.6421e-02   2.3925e-01   1.7468e-01  -4.7191e-02   1.9287e-01
    1.2043e+00   2.7622e-01   2.6765e-02   2.3074e-01   1.6551e-01  -4.7598e-02   1.8382e-01
    1.2217e+00   2.6803e-01   2.7087e-02   2.2274e-01   1.5686e-01  -4.7956e-02   1.7532e-01
    1.2392e+00   2.6032e-01   2.7401e-02   2.1522e-01   1.4872e-01  -4.8266e-02   1.6731e-01
    1.2566e+00   2.5305e-01   2.7694e-02   2.0814e-01   1.4103e-01  -4.8527e-02   1.5976e-01
    1.2741e+00   2.4618e-01   2.7961e-02   2.0147e-01   1.3377e-01  -4.8736e-02   1.5264e-01
    1.2915e+00   2.3969e-01   2.8204e-02   1.9518e-01   1.2692e-01  -4.8897e-02   1.4592e-01
    1.3090e+00   2.3355e-01   2.8425e-02   1.8925e-01   1.2044e-01  -4.9003e-02   1.3958e-01
    1.3265e+00   2.2774e-01   2.8609e-02   1.8364e-01   1.1431e-01  -4.9064e-02   1.3359e-01
    1.3439e+00   2.2223e-01   2.8757e-02   1.7835e-01   1.0851e-01  -4.9082e-02   1.2791e-01
    1.3614e+00   2.1700e-01   2.8874e-02   1.7334e-01   1.0299e-01  -4.9064e-02   1.2253e-01
    1.3788e+00   2.1204e-01   2.8952e-02   1.6861e-01   9.7767e-02  -4.9019e-02   1.1743e-01
    1.3963e+00   2.0733e-01   2.8987e-02   1.6412e-01   9.2813e-02  -4.8957e-02   1.1260e-01
    1.4137e+00   2.0286e-01   2.9001e-02   1.5989e-01   8.8102e-02  -4.8873e-02   1.0802e-01
    1.4312e+00   1.9862e-01   2.8987e-02   1.5589e-01   8.3627e-02  -4.8777e-02   1.0367e-01
    1.4486e+00   1.9460e-01   2.8949e-02   1.5210e-01   7.9381e-02  -4.8671e-02   9.9548e-02
    1.4661e+00   1.9079e-01   2.8895e-02   1.4853e-01   7.5339e-02  -4.8552e-02   9.5635e-02
    1.4835e+00   1.8718e-01   2.8820e-02   1.4514e-01   7.1497e-02  -4.8422e-02   9.1928e-02
    1.5010e+00   1.8374e-01   2.8719e-02   1.4193e-01   6.7839e-02  -4.8274e-02   8.8405e-02
    1.5184e+00   1.8048e-01   2.8602e-02   1.3888e-01   6.4347e-02  -4.8111e-02   8.5057e-02
    1.5359e+00   1.7737e-01   2.8457e-02   1.3597e-01   6.1014e-02  -4.7932e-02   8.1874e-02
    1.5533e+00   1.7440e-01   2.8288e-02   1.3319e-01   5.7824e-02  -4.7735e-02   7.8844e-02
    1.5708e+00   1.7157e-01   2.8101e-02   1.3053e-01   5.4760e-02  -4.7518e-02   7.5951e-02
    1.5882e+00   1.6887e-01   2.7894e-02   1.2798e-01   5.1819e-02  -4.7282e-02   7.3188e-02
    1.6057e+00   1.6629e-01   2.7671e-02   1.2554e-01   4.8986e-02  -4.7023e-02   7.0545e-02
    1.6232e+00   1.6382e-01   2.7428e-02   1.2319e-01   4.6251e-02  -4.6743e-02   6.8012e-02
    1.6406e+00   1.6145e-01   2.7169e-02   1.2093e-01   4.3609e-02  -4.6438e-02   6.5578e-02
    1.6581e+00   1.5916e-01   2.6884e-02   1.1874e-01   4.1054e-02  -4.6107e-02   6.3239e-02
    1.6755e+00   1.5696e-01   2.6586e-02   1.1662e-01   3.8574e-02  -4.5751e-02   6.0987e-02
    1.6930e+00   1.5484e-01   2.6270e-02   1.1457e-01   3.6168e-02  -4.5367e-02   5.8811e-02
    1.7104e+00   1.5279e-01   2.5938e-02   1.1259e-01   3.3820e-02  -4.4952e-02   5.6700e-02
    1.7279e+00   1.5081e-01   2.5585e-02   1.1067e-01   3.1537e-02  -4.4515e-02   5.4655e-02
    1.7453e+00   1.4889e-01   2.5213e-02   1.0880e-01   2.9312e-02  -4.4052e-02   5.2670e-02
    1.7628e+00   1.4703e-01   2.4826e-02   1.0699e-01   2.7133e-02  -4.3565e-02   5.0731e-02
    1.7802e+00   1.4523e-01   2.4425e-02   1.0523e-01   2.5006e-02  -4.3061e-02   4.8842e-02
    1.7977e+00   1.4349e-01   2.4010e-02   1.0353e-01   2.2921e-02  -4.2538e-02   4.6997e-02
    1.8151e+00   1.4181e-01   2.3580e-02   1.0188e-01   2.0887e-02  -4.2004e-02   4.5202e-02
    1.8326e+00   1.4018e-01   2.3137e-02   1.0028e-01   1.8892e-02  -4.1457e-02   4.3445e-02
    1.8500e+00   1.3861e-01   2.2686e-02   9.8732e-02   1.6935e-02  -4.0900e-02   4.1724e-02
    1.8675e+00   1.3709e-01   2.2220e-02   9.7224e-02   1.5021e-02  -4.0335e-02   4.0044e-02
    1.8850e+00   1.3563e-01   2.1739e-02   9.5767e-02   1.3155e-02  -3.9769e-02   3.8413e-02
    1.9024e+00   1.3423e-01   2.1247e-02   9.4360e-02   1.1322e-02  -3.9199e-02   3.6817e-02
    1.9199e+00   1.3288e-01   2.0743e-02   9.2995e-02   9.5260e-03  -3.8628e-02   3.5254e-02
    1.9373e+00   1.3158e-01   2.0230e-02   9.1670e-02   7.7634e-03  -3.8057e-02   3.3728e-02
    1.9548e+00   1.3033e-01   1.9714e-02   9.0392e-02   6.0233e-03  -3.7480e-02   3.2220e-02
    1.9722e+00   1.2912e-01   1.9187e-02   8.9146e-02   4.3187e-03  -3.6902e-02   3.0745e-02
    1.9897e+00   1.2795e-01   1.8653e-02   8.7934e-02   2.6456e-03  -3.6322e-02   2.9299e-02
    2.0071e+00   1.2683e-01   1.8113e-02   8.6767e-02   9.9866e-04  -3.5742e-02   2.7877e-02
    2.0246e+00   1.2574e-01   1.7567e-02   8.5633e-02  -6.2715e-04  -3.5158e-02   2.6468e-02
    2.0420e+00   1.2468e-01   1.7014e-02   8.4529e-02  -2.2282e-03  -3.4570e-02   2.5078e-02
    2.0595e+00   1.2364e-01   1.6452e-02   8.3447e-02  -3.7948e-03  -3.3981e-02   2.3718e-02
    2.0769e+00   1.2263e-01   1.5881e-02   8.2393e-02  -5.3324e-03  -3.3391e-02   2.2381e-02
    2.0944e+00   1.2166e-01   1.5306e-02   8.1375e-02  -6.8540e-03  -3.2802e-02   2.1061e-02
    2.1118e+00   1.2071e-01   1.4727e-02   8.0372e-02  -8.3507e-03  -3.2213e-02   1.9765e-02
    2.1293e+00   1.1978e-01   1.4152e-02   7.9383e-02  -9.8319e-03  -3.1612e-02   1.8487e-02
    2.1468e+00   1.1886e-01   1.3581e-02   7.8396e-02  -1.1300e-02  -3.0993e-02   1.7223e-02
    2.1642e+00   1.1795e-01   1.3010e-02   7.7412e-02  -1.2749e-02  -3.0354e-02   1.5980e-02
    2.1817e+00   1.1706e-01   1.2442e-02   7.6440e-02  -1.4187e-02  -2.9698e-02   1.4752e-02
    2.1991e+00   1.1618e-01   1.1876e-02   7.5473e-02  -1.5598e-02  -2.9029e-02   1.3547e-02
    2.2166e+00   1.1531e-01   1.1310e-02   7.4514e-02  -1.6983e-02  -2.8344e-02   1.2364e-02
    2.2340e+00   1.1447e-01   1.0750e-02   7.3584e-02  -1.8353e-02  -2.7646e-02   1.1190e-02
    2.2515e+00   1.1366e-01   1.0198e-02   7.2679e-02  -1.9711e-02  -2.6932e-02   1.0022e-02
    2.2689e+00   1.1287e-01   9.6475e-03   7.1790e-02  -2.1046e-02  -2.6203e-02   8.8736e-03
    2.2864e+00   1.1209e-01   9.0946e-03   7.0913e-02  -2.2362e-02  -2.5458e-02   7.7359e-03
    2.3038e+00   1.1136e-01   8.5377e-03   7.0078e-02  -2.3671e-02  -2.4714e-02   6.6054e-03
    2.3213e+00   1.1067e-01   7.9754e-03   6.9284e-02  -2.4965e-02  -2.3980e-02   5.4844e-03
    2.3387e+00   1.1005e-01   7.4176e-03   6.8549e-02  -2.6252e-02  -2.3267e-02   4.3698e-03
    2.3562e+00   1.0948e-01   6.8665e-03   6.7857e-02  -2.7515e-02  -2.2571e-02   3.2785e-03
    2.3736e+00   1.0897e-01   6.3289e-03   6.7219e-02  -2.8760e-02  -2.1894e-02   2.1991e-03
    2.3911e+00   1.0853e-01   5.8104e-03   6.6647e-02  -2.9996e-02  -2.1234e-02   1.1257e-03
    2.4086e+00   1.0813e-01   5.3062e-03   6.6118e-02  -3.1197e-02  -2.0585e-02   7.3005e-05
    2.4260e+00   1.0779e-01   4.8129e-03   6.5653e-02  -3.2369e-02  -1.9954e-02  -9.6092e-04
    2.4435e+00   1.0753e-01   4.3332e-03   6.5271e-02  -3.3530e-02  -1.9352e-02  -1.9922e-03
    2.4609e+00   1.0734e-01   3.8692e-03   6.4963e-02  -3.4675e-02  -1.8777e-02  -3.0179e-03
    2.4784e+00   1.0720e-01   3.4160e-03   6.4711e-02  -3.5790e-02  -1.8227e-02  -4.0230e-03
    2.4958e+00   1.0715e-01   2.9689e-03   6.4542e-02  -3.6895e-02  -1.7712e-02  -5.0170e-03
    2.5133e+00   1.0718e-01   2.5201e-03   6.4442e-02  -3.7985e-02  -1.7236e-02  -5.9921e-03
    2.5307e+00   1.0729e-01   2.0731e-03   6.4407e-02  -3.9069e-02  -1.6797e-02  -6.9554e-03
    2.5482e+00   1.0749e-01   1.6347e-03   6.4433e-02  -4.0166e-02  -1.6383e-02  -7.9108e-03
    2.5656e+00   1.0776e-01   1.2069e-03   6.4498e-02  -4.1281e-02  -1.5972e-02  -8.8628e-03
    2.5831e+00   1.0809e-01   7.8616e-04   6.4588e-02  -4.2399e-02  -1.5569e-02  -9.7894e-03
    2.6005e+00   1.0850e-01   3.7092e-04   6.4721e-02  -4.3522e-02  -1.5186e-02  -1.0689e-02
    2.6180e+00   1.0904e-01  -3.9518e-05   6.4951e-02  -4.4673e-02  -1.4841e-02  -1.1584e-02
    2.6354e+00   1.0971e-01  -4.4440e-04   6.5287e-02  -4.5830e-02  -1.4552e-02  -1.2463e-02
    2.6529e+00   1.1051e-01  -8.3137e-04   6.5729e-02  -4.6991e-02  -1.4306e-02  -1.3323e-02
    2.6704e+00   1.1140e-01  -1.1921e-03   6.6243e-02  -4.8145e-02  -1.4073e-02  -1.4156e-02
    2.6878e+00   1.1237e-01  -1.5295e-03   6.6823e-02  -4.9297e-02  -1.3836e-02  -1.4972e-02
    2.7053e+00   1.1342e-01  -1.8551e-03   6.7478e-02  -5.0444e-02  -1.3595e-02  -1.5777e-02
    2.7227e+00   1.1455e-01  -2.1888e-03   6.8225e-02  -5.1561e-02  -1.3367e-02  -1.6556e-02
    2.7402e+00   1.1579e-01  -2.5433e-03   6.9114e-02  -5.2644e-02  -1.3173e-02  -1.7331e-02
    2.7576e+00   1.1711e-01  -2.9123e-03   7.0141e-02  -5.3672e-02  -1.3004e-02  -1.8103e-02
    2.7751e+00   1.1844e-01  -3.2693e-03   7.1246e-02  -5.4627e-02  -1.2803e-02  -1.8871e-02
    2.7925e+00   1.1971e-01  -3.5882e-03   7.2356e-02  -5.5508e-02  -1.2495e-02  -1.9641e-02
    2.8100e+00   1.2085e-01  -3.8718e-03   7.3403e-02  -5.6321e-02  -1.2014e-02  -2.0416e-02
    2.8274e+00   1.2189e-01  -4.1499e-03   7.4394e-02  -5.7088e-02  -1.1352e-02  -2.1203e-02
    2.8449e+00   1.2284e-01  -4.4415e-03   7.5327e-02  -5.7795e-02  -1.0515e-02  -2.1975e-02
    2.8623e+00   1.2371e-01  -4.7387e-03   7.6202e-02  -5.8442e-02  -9.5014e-03  -2.2723e-02
    2.8798e+00   1.2457e-01  -5.0389e-03   7.7090e-02  -5.9035e-02  -8.3235e-03  -2.3450e-02
    2.8972e+00   1.2548e-01  -5.3432e-03   7.8069e-02  -5.9549e-02  -6.9938e-03  -2.4146e-02
    2.9147e+00   1.2639e-01  -5.6313e-03   7.9142e-02  -5.9977e-02  -5.4724e-03  -2.4852e-02
    2.9322e+00   1.2730e-01  -5.8674e-03   8.0340e-02  -6.0387e-02  -3.6884e-03  -2.5659e-02
    2.9496e+00   1.2815e-01  -5.9672e-03   8.1610e-02  -6.0866e-02  -1.5151e-03  -2.6646e-02
    2.9671e+00   1.2890e-01  -5.7872e-03   8.2800e-02  -6.1626e-02   1.2444e-03  -2.7882e-02
    2.9845e+00   1.2969e-01  -5.1890e-03   8.3768e-02  -6.3012e-02   4.6936e-03  -2.9373e-02
    3.0020e+00   1.3113e-01  -4.0994e-03   8.4634e-02  -6.5532e-02   8.6071e-03  -3.1057e-02
    3.0194e+00   1.3432e-01  -2.5591e-03   8.5871e-02  -6.9770e-02   1.2305e-02  -3.2790e-02
    3.0369e+00   1.4061e-01  -8.0056e-04   8.8237e-02  -7.6168e-02   1.4737e-02  -3.4345e-02
    3.0543e+00   1.5110e-01   7.7478e-04   9.2520e-02  -8.4858e-02   1.4891e-02  -3.5524e-02
    3.0718e+00   1.6578e-01   1.6724e-03   9.8957e-02  -9.5259e-02   1.2267e-02  -3.6089e-02
    3.0892e+00   1.8278e-01   1.4960e-03   1.0687e-01  -1.0577e-01   7.3653e-03  -3.5867e-02
    3.1067e+00   1.9736e-01   5.5354e-04   1.1398e-01  -1.1379e-01   2.2487e-03  -3.5108e-02
    3.1241e+00   2.0333e-01  -1.6825e-07   1.1696e-01  -1.1684e-01   5.7126e-06  -3.4674e-02
    3.1416e+00   2.0980e-01            0   1.2006e-01  -1.2003e-01            0  -3.4418e-02];
    Pmat(:,:,2)=Pmat(:,:,1);
    [thetak,wthetak]=lgwt(1500,0,pi);
    for n=1:length(SSA_p_TYPE)
        THETA=Pmat(:,1,n);
        % check whether the integral result is euqal to 2 and normalize anyway        
        P11_nmt=interp1(Pmat(:,1),Pmat(:,2,n),thetak,'spline');
        sumP11=P11_nmt.*sin(thetak)*wthetak'    
        Pmat(:,2:end,n)=Pmat(:,2:end,n)*(2/sumP11);
    
        F11_nm=Pmat(:,2,n);     
        F21_nm=Pmat(:,3,n);         
        F22_nm=Pmat(:,4,n); 
        F33_nm=Pmat(:,5,n);             
        F43_nm=Pmat(:,6,n);             
        F44_nm=Pmat(:,7,n); 

        %plot(THETA,log10(F21_nm))

        % convert theta points to 
        P11_nm_TYPE(n,:)=interp1(THETA,F11_nm,theta,'spline');
        P21_nm_TYPE(n,:)=interp1(THETA,F21_nm,theta,'spline');
        P22_nm_TYPE(n,:)=interp1(THETA,F22_nm,theta,'spline');
        P33_nm_TYPE(n,:)=interp1(THETA,F33_nm,theta,'spline');
        P43_nm_TYPE(n,:)=interp1(THETA,F43_nm,theta,'spline');
        P44_nm_TYPE(n,:)=interp1(THETA,F44_nm,theta,'spline');
       
    end
end

% Delta truncation
if delta_use==1
    for n=1:Aerosol_Type
        if SSA_p_TYPE(n)~=0
            %[OD(:,n+2),SSA_p_TYPE(n),P11_nm_TYPE(n,:),P21_nm_TYPE(n,:),P22_nm_TYPE(n,:),P33_nm_TYPE(n,:),P43_nm_TYPE(n,:),P44_nm_TYPE(n,:)]=Delta_truncation_PA(OD(:,n+2),SSA_p_TYPE(n),theta,P11_nm_TYPE(n,:),P21_nm_TYPE(n,:),P22_nm_TYPE(n,:),P33_nm_TYPE(n,:),P43_nm_TYPE(n,:),P44_nm_TYPE(n,:),trunc_angle);
            [OD(:,n+2),SSA_p_TYPE(n),P11_nm_TYPE(n,:),P21_nm_TYPE(n,:),P22_nm_TYPE(n,:),P33_nm_TYPE(n,:),P43_nm_TYPE(n,:),P44_nm_TYPE(n,:)]=Delta_truncation_PA_Imp(OD(:,n+2),SSA_p_TYPE(n),theta,P11_nm_TYPE(n,:),P21_nm_TYPE(n,:),P22_nm_TYPE(n,:),P33_nm_TYPE(n,:),P43_nm_TYPE(n,:),P44_nm_TYPE(n,:));
        end
    %         sumP11=0;
    %         n=2
    %         for thetan=1:length(theta)-1
    %             sumP11=sumP11+P11_nm_TYPE(n,thetan)*[theta(thetan+1)-theta(thetan)]*sin(theta(thetan));
    %         end
    %         % check whether the integral result is euqal to 2 (apprixmately)
    %         sumP11                 
        
    end
end
Lay_num=OD(end,1);
SSA_p(Lay_num)=0;   
AF_p(Lay_num,Aerosol_Type)=0;
for L=1:Lay_num
    %  and phase matrix for mixed aerosol scattering
    ODA(L)=sum(OD(L,3:end));
    if ODA(L)~=0
        %% mixed aerosol scattering albedo
        scaOD=sum(OD(L,3:end).*SSA_p_TYPE);
        if scaOD~=0
            SSA_p(L)=sum(OD(L,3:end).*SSA_p_TYPE)./ODA(L);
            %% fraction of differnet types of aerosol scattering 
            AF_p(L,:)=OD(L,3:end).*SSA_p_TYPE/[SSA_p(L)*ODA(L)];
        end
    end        
end

L=1;
MiemixtypeforL(L)=L;
Miemixtype(L)=L;    
Aerosoleq0=zeros(1,Aerosol_Type);
for L=2:Lay_num
    if (isequal(AF_p(L,:),Aerosoleq0)==1)|(abs(AF_p(L,:)-AF_p(L-1,:))<1e-8);
        MiemixtypeforL(L)=MiemixtypeforL(L-1);
    else
        MiemixtypeforL(L)=MiemixtypeforL(L-1)+1;
        Miemixtype(MiemixtypeforL(L-1)+1)=L;
    end       
end   

SSA_m(1:Lay_num)=1;   
OD_mixed=[OD(:,1:2),ODA'];
% optical thickness of each layer counted from the TOA
TAUALL(1)=OD_mixed(1,2)+OD_mixed(1,3);
for L=2:OD_mixed(end,1)
    TAUALL(L)=TAUALL(L-1)+[OD_mixed(L,2)+OD_mixed(L,3)];
end
% fraction of Rayleigh scattering
RAYF=OD_mixed(:,2)'.*SSA_m./[OD_mixed(:,2)'.*SSA_m+OD_mixed(:,3)'.*SSA_p];
% Single scattering albedo
OMEGA=[OD_mixed(:,2)'.*SSA_m+OD_mixed(:,3)'.*SSA_p]./(OD_mixed(:,2)+OD_mixed(:,3))'; 
    
LMiemixtype=length(Miemixtype);



%% Specification of input parameters for computation 
deta=0.0; % depolarization factor: 0 for isotropic Rayleigh scattering (2.908633816103590d-2 for standard Earth atmosphere)

%% RPV surface model parameters for non-polarizing reflection part
RPV_b=-0.5; % RPV model parameter: typically between [-1,1] 
RPV_k=1.5;  % RPV model parameter: typically between [0,2]
RPV_r=0.0;    %  RPV model parameter: typically between [0,infinity] (set this value to be zero for black surface) ("0.17096" used in F.Xu's paper)

%% Cox-munk's model parameters for polarizing reflection part  
nn=1.334;   % refrative index of the polarizing surface
vv=2.0;       % surface wind speed (m/s)
epsirol=0.0;% Weight of polarization contribution; % epsirol is typically between [0,1] (also set this value to be zero for black surface)


M_Max=31;% Maximum Fourier series expansion number
AngGaussNum=16;% % Stream number
[XMU,WTMU]=lgwt(AngGaussNum,0,1); 
%% viewing azimuthal angles
faipfai0=[0,pi];%[0,3*pi/4,pi]; 
Lfaipfai0=length(faipfai0);
%% Extra incidence angle in addition to the Gauss angular points
XMU0_ADD=cos([10,60]*pi/180);%if no extra incidence angle specified use "XMU0_ADD=[]" here.
%% Extra viewing angle in addition to the  Gauss angular points 
XMUe_ADD=cos([0:1:90]*pi/180); %if no extra viewing angle specified use "XMUe_ADD=[]" here.
XMUe_ADD(find(XMUe_ADD==1))=1-1.0E-10;
dex=0;
for n=1:length(XMUe_ADD)
    dex=find(abs(XMU-XMUe_ADD(n))<1.0e-9);
    if dex~=0
        display(['The ',num2str(n),'-th',' viewing angle in "XMUe_ADD" is already contained in "XMU-array" and will be computed, please delete it.'])
        return
    end
end
dex=0;
for n=1:length(XMU0_ADD)
    dex=find(abs(XMU-XMU0_ADD(n))<1.0e-9);
    if dex~=0
        display(['The ',num2str(n),'-th',' incident angle in "XMU0_ADD" is already contained in "XMU-array" and will be computed, please delete it.'])
        return
    end
end
XMU0ALL=[XMU,XMU0_ADD];
XMUeALL=[XMU,XMUe_ADD];

NUMLAY_GROUP=3;    % Number of sublayers in each group (should be 2<=NUMLAY_GROUP <=Lay_num) 
NMultiplication=3; % Number of expansion terms in approximating inverse of "I-Q", where I is identity matrix
L_max=1313;        % Moment number for FFT decomposition of phase matrix
xi_number=1500;    % Gauss point number of -1<=ksi<=1 for alfa, beta, gamma and delta integral (refer to Siewert's Paper 1982)
%% =========  Kernal Markov chain computation part ========== %%
Computation_mode='vector2'; % 'vector1' for the full vector radiative transfer code with {I, Q, U, V} computation
                            % 'vector2' for the quasi-vector RT code with {I, Q, U} computation 
                            
if isequal(Computation_mode,'vector1')==1
    DStokes=4;    
    [IS,IM]=VectorMarkovchainStokes34(deta,RPV_b,RPV_k,RPV_r,nn,vv,epsirol,M_Max,XMU,WTMU,XMU0ALL,XMUeALL,faipfai0,Lfaipfai0,NUMLAY_GROUP,NMultiplication,L_max,xi_number,DStokes);
    [IS(1:8,7,1),IM(1:8,7,1)]
    QU_cal=1;
    V_cal=1;
end
if isequal(Computation_mode,'vector2')==1
    DStokes=3;    
    [IS,IM]=VectorMarkovchainStokes34(deta,RPV_b,RPV_k,RPV_r,nn,vv,epsirol,M_Max,XMU,WTMU,XMU0ALL,XMUeALL,faipfai0,Lfaipfai0,NUMLAY_GROUP,NMultiplication,L_max,xi_number,DStokes);
    [IS(1:8,7,1),IM(1:8,7,1)]
    
    QU_cal=1;
    V_cal=0;
end

% figure(1)
% axes1 = axes('Fontsize',14,'box','on') 
%LE=(length(XMU)+1)/2;
LE=length(XMU0ALL);
[XMUeALL2,k]=sort(XMUeALL);

if isequal(Computation_mode,'vector1')==1
    ISMview=[];
    QSMview=[];
    USMview=[];
    VSMview=[];    
    for jjj=1:Lfaipfai0  
        % single scattered light
        Isingle=IS(1:4:end,LE,jjj);
        Qsingle=IS(2:4:end,LE,jjj);
        Usingle=IS(3:4:end,LE,jjj);
        Vsingle=IS(4:4:end,LE,jjj);

        % multiple scattered light
        Imultiple=IM(1:4:end,LE,jjj);
        Qmultiple=IM(2:4:end,LE,jjj);
        Umultiple=IM(3:4:end,LE,jjj); 
        Vmultiple=IM(4:4:end,LE,jjj);

        for n=1:length(k)
            ISM2(n)=Isingle(k(n))+Imultiple(k(n));
            QSM2(n)=Qsingle(k(n))+Qmultiple(k(n));
            USM2(n)=Usingle(k(n))+Umultiple(k(n));
            VSM2(n)=Vsingle(k(n))+Vmultiple(k(n));        
        end
        XMUview=XMUeALL2;  
        LXMUview=length(XMUview);        
        ISMview(jjj,1:LXMUview)=ISM2;
        QSMview(jjj,1:LXMUview)=QSM2;
        USMview(jjj,1:LXMUview)=USM2;
        VSMview(jjj,1:LXMUview)=VSM2; 
    end    
end

if isequal(Computation_mode,'vector2')==1
    ISMview=[];
    QSMview=[];
    USMview=[];
    VSMview=[];
    for jjj=1:Lfaipfai0  
        % single scattered light
        Isingle=IS(1:3:end,LE,jjj);
        Qsingle=IS(2:3:end,LE,jjj);
        Usingle=IS(3:3:end,LE,jjj);

        % multiple scattered light
        Imultiple=IM(1:3:end,LE,jjj);
        Qmultiple=IM(2:3:end,LE,jjj);
        Umultiple=IM(3:3:end,LE,jjj); 

        for n=1:length(k)
            ISM2(n)=Isingle(k(n))+Imultiple(k(n));
            QSM2(n)=Qsingle(k(n))+Qmultiple(k(n));
            USM2(n)=Usingle(k(n))+Umultiple(k(n));       
        end


        XMUview=XMUeALL2;  
        LXMUview=length(XMUview);        
        ISMview(jjj,1:LXMUview)=ISM2;
        QSMview(jjj,1:LXMUview)=QSM2;
        USMview(jjj,1:LXMUview)=USM2;

    end
end

Nstep=10;
%figure('position',[300,300,650,500])
%axes1 = axes('Fontsize',14,'box','on') 
%plot([-acos(XMUview),fliplr(acos(XMUview))]*180/pi,[ISMview(2,:),fliplr(ISMview(1,:))],'linestyle','-','color','b','linewidth',1.0);
%hold on
%plot([-acos(XMUview),fliplr(acos(XMUview))]*180/pi,-10*[QSMview(2,:),fliplr(QSMview(1,:))],'linestyle','-','color','r','linewidth',1.0);
%hold on
%legend('Markov Chain: I ','Markov Chain: -10Q')
output_array = transpose([[-acos(XMUview),fliplr(acos(XMUview))]*180/pi;[ISMview(2,:),fliplr(ISMview(1,:))];[QSMview(2,:),fliplr(QSMview(1,:))];[USMview(2,:),fliplr(USMview(1,:))]])
output_file=parameters.output_file
csvwrite(output_file, output_array)

%xlim([-90,90])
%xTickLabels=[-90:30:90];
%set(gca,'Xtick', xTickLabels);
%xlabel('\Theta (deg)','Fontsize',14,'Fontname','Times New Roman')
%ylabel('\pi*S/\mu_{0} ','Fontsize',14,'Fontname','Times New Roman')
%% ========================================================== %% 
