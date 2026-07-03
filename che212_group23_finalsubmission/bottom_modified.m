clc;
clear;
close all;
%% Given Data
V=36;%in volts
I=0.22;%in ampere
v=[7.2,11.3,12.6];%in m/s
Ts=[40.46+273.15,39.16+273.15,38.7+273.15];%in K
T_inf=[30.2+273.15,30.9+273.15,30.9+273.15];%in K
d=0.03986;L=0.2;%in metres
Pr=0.71;%Prandtl number
K=0.027;%in W/mK
neu=1.7*10^-5;
sigma=5.67*10^-8;eps=0.85;
g=9.8;%in m2/s
%% Experimental heat tranfer cofficient
for i=1:length(Ts)
    h_exp(i)=(V*I)./((3.14*d*L).*(Ts(i)-T_inf(i)));
end
%% Forced heat tranfer cofficient 
for j=1:length(Ts)
    Re(j)=(v(j).*L)./neu;%Reynold's Number
    Nu(j)=0.664.*(Pr.^(1/3)).*(Re(j).^(1/2));%Nusselt correlation
    hf(j)=(Nu(j).*K)./L;
end
%% Radiation heat tranfer cofficient
for j=1:length(Ts)
    hr(j)=sigma.*eps.*(Ts(j)+T_inf(j)).*((Ts(j).^2)+T_inf(j).^2);
end
%% Natural heat tranfer cofficient
for j=1:length(Ts)
    Tf(j)=(Ts(j)+T_inf(j))./2;
    beta(j)=1./Tf(j);
    Gr(j)=(g.*(Ts(j)-T_inf(j)).*(0.008))./((2.89*10^-10).*Tf(j));
    Nu_n(j)=0.59.*(Gr(j).*Pr).^(1/4);
    hn(j)=(Nu_n(j).*K)./L;
end
%% Overall heat tranfer cofficient & error measurement
for j=1:length(Ts)
    Nu_mix(j)=(Nu_n(j).^3+Nu(j).^3).^(1/3);
    h_conv(j)=(Nu_mix(j).*K)./L;
    h_theo(j)=hr(j)+h_conv(j);
    h_t(j)=hf(j)+hr(j)+hn(j);
    error(j)=(1-(h_exp(j)./h_theo(j))).*100;
end
