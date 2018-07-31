function [OD,SSA_p,P11_nm,P21_nm,P22_nm,P33_nm,P43_nm,P44_nm]=Delta_truncation_PA(OD,SSA_p,theta,P11_nm,P21_nm,P22_nm,P33_nm,P43_nm,P44_nm,trunc_angle)
%determine the cutting-angle by fixing an angle value 
for n=1:length(theta)
    if theta(n)>trunc_angle*pi/180;
        delta_N=n;
        kn2=[log10(P11_nm(n))-log10(P11_nm(n-1))]/[theta(n)-theta(n-1)]; 
        delta_k=kn2;
        delta_need=1;
        break
    end
end
% applying the delta-approximation 
if delta_need==1
    theta_d=theta(1:delta_N);
    theta_N=theta(delta_N);
    % a: replacing P11 in [1, theta(delta_N)] by parabolic f=a*x^2+b
    a=delta_k/2/theta_N; 
    b=log10(P11_nm(delta_N))-a*theta_N^2;
    P11_nm_d=10.^(a*theta_d.^2+b); 

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

    OD=(1-trunc_ratio*SSA_p)*OD;
    SSA_p=[(1-trunc_ratio)*SSA_p]/[1-trunc_ratio*SSA_p]; %NAKAJIMA and ASANO's Eq.(1)


%         plot(theta,log10(P11_nm_ori))
%         hold on
%         plot(theta,log10(P11_nm),'r') 

    FP=P11_nm./P11_nm_ori;
    P22_nm=P22_nm./P11_nm_ori.*P11_nm;
    P21_nm=P21_nm./P11_nm_ori.*P11_nm;
    P33_nm=P33_nm./P11_nm_ori.*P11_nm;
    P43_nm=P43_nm./P11_nm_ori.*P11_nm;
    P44_nm=P44_nm./P11_nm_ori.*P11_nm;

end

