function [OD,SSA_p,P11_nm,P21_nm,P22_nm,P33_nm,P43_nm,P44_nm,f_trunc]=Delta_truncation_PA_Imp(OD,SSA_p,theta,P11_nm,P21_nm,P22_nm,P33_nm,P43_nm,P44_nm)


dtheta=1*pi/180; %% angular stepsize to get the dP/dtheta 
P11value_cut=10;
kkk=find(P11_nm<P11value_cut);

P11_trunval=2*P11_nm(kkk(1));
delta_N=find(P11_nm<P11_trunval);
delta_N=delta_N(1);
theta_trun=theta(delta_N);
theta2=theta_trun-dtheta;
f_trunc=0;
if theta2>0
    delta_N2=find(theta>=theta2);
    delta_N2=delta_N2(1);
    log10use=1;
    if log10use==1 
        delta_k=[log10(P11_nm(delta_N))-log10(P11_nm(delta_N2))]/[theta(delta_N)-theta(delta_N2)];   % my strategy
    else
        delta_k=[P11_nm(delta_N)-P11_nm(delta_N2)]/[theta(delta_N)-theta(delta_N2)];  % John's strategy
    end
    theta_d=theta(1:delta_N);
    theta_N=theta(delta_N);
    % a: replacing P11 in [1, theta(delta_N)] by parabolic f=a*x^2+b
    a=delta_k/2/theta_N;
    if log10use==1 
        b=log10(P11_nm(delta_N))-a*theta_N^2;
        P11_nm_d=10.^(a*theta_d.^2+b); 
    else
        b=P11_nm(delta_N)-a*theta_N^2;
        P11_nm_d=a*theta_d.^2+b; 
    end

    % b: correcting the aerosol optical depth when delta-approximation is to be applied
    [x,wx]=lgwt(50,0,theta_N);
    % original
    y=interp1(theta,P11_nm,x,'spline');
    IntP11=sum(wx.*y.*sin(x));   
    % approximated
    y=interp1(theta_d,P11_nm_d,x,'spline');
    IntP11_d=sum(wx.*y.*sin(x));


    AA=IntP11-IntP11_d;
    trunc_ratio=AA/2; % the truncation ratio

    P11_nm_ori=P11_nm;
    P11_nm(1:delta_N)=P11_nm_d;
    P11_nm=P11_nm/(1-trunc_ratio); 

    f_trunc=trunc_ratio*SSA_p;

    OD=(1-trunc_ratio*SSA_p)*OD; %NAKAJIMA and ASANO's Eq.(2)        
    SSA_p=[(1-trunc_ratio)*SSA_p]/[1-trunc_ratio*SSA_p]; %NAKAJIMA and ASANO's Eq.(1)


    % check whether the integral result is euqal to 2
    [thetak,wthetak]=lgwt(1500,0,pi);
    P11_nmt=interp1(theta,P11_nm,thetak,'spline');
    sumP11=P11_nmt.*sin(thetak)*wthetak';

    [thetak,wthetak]=lgwt(1500,0,pi);
    P11_nmt=interp1(theta,P11_nm_ori,thetak,'spline');
    sumP11_ori=P11_nmt.*sin(thetak)*wthetak';

    FP=P11_nm./P11_nm_ori;
    P22_nm=P22_nm.*FP;
    P21_nm=P21_nm.*FP;
    P33_nm=P33_nm.*FP;
    P43_nm=P43_nm.*FP;
    P44_nm=P44_nm.*FP;
    %     plot(theta*180/pi,log10(P11_nm_ori))
    %     hold on
    %     plot(theta*180/pi,log10(P11_nm),'r')     
    %     a=1  
end

      
        