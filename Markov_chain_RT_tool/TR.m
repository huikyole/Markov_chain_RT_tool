function [Tm_bar,Rm_bar,Pnm_bar]=TR(miu,L_max,M_Max)
Lmiu=length(miu);
Pnm_bar(1:L_max+1,1:Lmiu,1:M_Max+1)=0;
Pnm(1:L_max+1,1:Lmiu)=0;
PI_mn(1:L_max+1,1:Lmiu)=0;
% tic
% Pnm=Pu_cal(miu,L_max); % Note: Pnm matrix structure: row--m;column--cosTheta;plane--n
% for m=0:M_Max
%     for L=m:L_max
%         Pnm_bar(L+1,:,m+1)=Pnm(m+1,:,L+1)/sqrt(prodc(L-m,L+m));
%     end
% end
% toc

% here to be added for Olga (3)
%==== a substitution to save time (referring to "PnmDPnm" function) ====
%tic
% Part I: m=0
m=0;
%% n=0
P_0=1;
Pnm(1,1:Lmiu)=P_0;
%% n=1
P_1=miu;
Pnm(2,1:Lmiu)=P_1;

P_0p(1:Lmiu)=0;

for n=1:L_max-1 
    P_1p=n*P_0+miu.*P_0p;  
    P_2=((2*n+1)/(n+1))*miu.*P_1-(n/(n+1))*P_0;
    P_0=P_1;P_1=P_2;P_0p=P_1p;
    Pnm(n+2,:)=P_2;
end
%Pnm=P_0n; %Attention that the values of PI_0n is designed from Eq.13,29-30, according to J.A.Lock,JOSA.A.1993,P696.
for L=m:L_max
    Pnm_bar(L+1,:,m+1)=Pnm(L+1,:)/sqrt(prodc(L-m,L+m));   
end

% Part II: m>=1, results corresponding to n=m..n_max 
for m=1:M_Max
    Factorialx=1;
    jj=[1:1:m];
    jj2=(2*jj-1);
    Factorialx=prod(jj2);

    PI_a(1:Lmiu)=0;
    PI_b(1:Lmiu)=Factorialx*[(sqrt(1-miu.^2)).^(m-1)];
    PI_mn(1,1:Lmiu)=PI_b;
    for n=m+1:L_max+1
       PI_mn(n-m+1,:)=[(2*(n-1)+1)/((n-1)+1-m)]*miu.*PI_b-[((n-1)+m)/((n-1)+1-m)]*PI_a;
       PI_a=PI_b;PI_b=PI_mn(n-m+1,:);
    end  

    %PI_mn=PI_store;
    %Following both effective in effect and important in conception 
    %PI_mn=(-1)^(m)*PI_mn; 
    %TAU_mn=(-1)^(m)*TAU_mn;
    for n=m:L_max+1
        Pnm(n-m+1,:)=PI_mn(n-m+1,:).*sqrt(1-miu.^2);
        Pnm_bar(n+1,:,m+1)=Pnm(n-m+1,:)/sqrt(prodc(n-m,n+m));
    end
end
%toc
%=====================================


% a) m=0
m=0;   
L=m;
Rm_bar(L+1,1:Lmiu,m+1)=0;  
Tm_bar(L+1,1:Lmiu,m+1)=0; 
L=m+1;
Rm_bar(L+1,1:Lmiu,m+1)=0;  
Tm_bar(L+1,1:Lmiu,m+1)=0; 
L=m+2;
Rm_bar(L+1,1:Lmiu,m+1)=sqrt(6)/4*(1-miu.^2);  %Min is wrong here in sign in his 2004 paper (Eq.21b) but correct in his 2009 paper (Eq.A15).
Tm_bar(L+1,1:Lmiu,m+1)=0; 
for L=2:L_max-1 % then calculate for L+1
    tmL=sqrt(((L+1)^2-m^2)*((L+1)^2-4))/(L+1);    
%     if m==L    %never happens           
%         Rm_bar(m+1,:,L+1+1)=(2*L+1)*(miu.*Rm_bar(m+1,:,L+1)-2*m/L/(L+1)*Tm_bar(m+1,:,L+1))/tmL;  
%         Tm_bar(m+1,:,L+1+1)=(2*L+1)*(miu.*Tm_bar(m+1,:,L+1)-2*m/L/(L+1)*Rm_bar(m+1,:,L+1))/tmL; 
%     else
        Rm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Rm_bar(L+1,:,m+1)-2*m/L/(L+1)*Tm_bar(L+1,:,m+1))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*Rm_bar(L-1+1,:,m+1)/tmL;  
        Tm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Tm_bar(L+1,:,m+1)-2*m/L/(L+1)*Rm_bar(L+1,:,m+1))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*Tm_bar(L-1+1,:,m+1)/tmL;            
%     end  
end    
% b) m=1
m=1;
L=m;
Rm_bar(L+1,1:Lmiu,m+1)=0;  
Tm_bar(L+1,1:Lmiu,m+1)=0; 
L=m+1;
Rm_bar(L+1,1:Lmiu,m+1)=-miu/2.*sqrt(1-miu.^2);
Tm_bar(L+1,1:Lmiu,m+1)=-1/2*sqrt(1-miu.^2);
for L=2:L_max-1 % then calculate for L+1
    tmL=sqrt(((L+1)^2-m^2)*((L+1)^2-4))/(L+1);
%     if m==L    %never happens           
%         Rm_bar(m+1,:,L+1+1)=(2*L+1)*(miu.*Rm_bar(m+1,:,L+1)-2*m/L/(L+1)*Tm_bar(m+1,:,L+1))/tmL;  
%         Tm_bar(m+1,:,L+1+1)=(2*L+1)*(miu.*Tm_bar(m+1,:,L+1)-2*m/L/(L+1)*Rm_bar(m+1,:,L+1))/tmL; 
%     else
        Rm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Rm_bar(L+1,:,m+1)-2*m/L/(L+1)*Tm_bar(L+1,:,m+1))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*Rm_bar(L-1+1,:,m+1)/tmL;  
        Tm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Tm_bar(L+1,:,m+1)-2*m/L/(L+1)*Rm_bar(L+1,:,m+1))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*Tm_bar(L-1+1,:,m+1)/tmL;            
%     end    
end
% % c) m>=2
% for m=2:M_Max
%     L=m;
%     Rm_bar(L+1,:,m+1)=sqrt(m*(m-1)/(m+1)/(m+2))*(1+miu.^2)./(1-miu.^2).*Pnm_bar(L+1,:,m+1);  
%     Tm_bar(L+1,:,m+1)=sqrt(m*(m-1)/(m+1)/(m+2))*(2*miu)./(1-miu.^2).*Pnm_bar(L+1,:,m+1);
%     for L=m:L_max-1 % then calculate for L+1
%         tmL=sqrt(((L+1)^2-m^2)*((L+1)^2-4))/(L+1);
%         if m==L            
%             Rm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Rm_bar(L+1,:,m+1)-2*m/L/(L+1)*Tm_bar(L+1,:,m+1))/tmL;  
%             Tm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Tm_bar(L+1,:,m+1)-2*m/L/(L+1)*Rm_bar(L+1,:,m+1))/tmL; 
%         else
%             Rm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Rm_bar(L+1,:,m+1)-2*m/L/(L+1)*Tm_bar(L+1,:,m+1))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*Rm_bar(L-1+1,:,m+1)/tmL;  
%             Tm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Tm_bar(L+1,:,m+1)-2*m/L/(L+1)*Rm_bar(L+1,:,m+1))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*Tm_bar(L-1+1,:,m+1)/tmL;            
%         end 
%     end
% end

% Replacing the above statements by removing if-statement
for m=2:M_Max
    L=m;
    Rm_bar(L+1,:,m+1)=sqrt(m*(m-1)/(m+1)/(m+2))*(1+miu.^2)./(1-miu.^2).*Pnm_bar(L+1,:,m+1);  
    Tm_bar(L+1,:,m+1)=sqrt(m*(m-1)/(m+1)/(m+2))*(2*miu)./(1-miu.^2).*Pnm_bar(L+1,:,m+1);

    tmL=sqrt(((L+1)^2-m^2)*((L+1)^2-4))/(L+1);    
    Rm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Rm_bar(L+1,:,m+1)-2*m/L/(L+1)*Tm_bar(L+1,:,m+1))/tmL;  
    Tm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Tm_bar(L+1,:,m+1)-2*m/L/(L+1)*Rm_bar(L+1,:,m+1))/tmL; 
            
    for L=m+1:L_max-1 % then calculate for L+1
        tmL=sqrt(((L+1)^2-m^2)*((L+1)^2-4))/(L+1);
        Rm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Rm_bar(L+1,:,m+1)-2*m/L/(L+1)*Tm_bar(L+1,:,m+1))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*Rm_bar(L-1+1,:,m+1)/tmL;  
        Tm_bar(L+1+1,:,m+1)=(2*L+1)*(miu.*Tm_bar(L+1,:,m+1)-2*m/L/(L+1)*Rm_bar(L+1,:,m+1))/tmL-sqrt((L^2-m^2)*(L^2-4))/L*Tm_bar(L-1+1,:,m+1)/tmL;            
    end
end

return
