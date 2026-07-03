clc;
clear;
close all;
%% Given Data
V=38;%in volts
I=0.22;%in ampere
v=[6.1,8.4,10.1];%in m/s
Ts=[39.43+273.15,37.77+273.15,37.27+273.15];%in K
T_inf=[31.3+273.15,31+273.15,31+273.15];%in K
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
    Re(j)=(v(j).*d)./neu;%Reynold's Number
    Nu(j)=0.193.*(Pr.^(1/3)).*(Re(j).^(0.618));%Hilpert correlation
    hf(j)=(Nu(j).*K)./d;
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
    hn(j)=(Nu_n(j).*K)./d;
end
%% Overall heat tranfer cofficient & error measurement
for j=1:length(Ts)
    Nu_mix(j)=(Nu_n(j).^3+Nu(j).^3).^(1/3);
    h_conv(j)=(Nu_mix(j).*K)./d;
    h_theo(j)=hr(j)+h_conv(j);
    error(j)=(1-(h_exp(j)./h_theo(j))).*100;
    end
