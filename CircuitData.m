close all;clear;
%circuit analysis of the piX excitor from Zurich Med Tech
k=.5; %coefficient of coupling
w0=2*pi*128e6;
L1=38.64e-9; %inductance values
L2=60.02e-9;
M=k*sqrt(L2*L1); %mutual inductance
C=1/(w0^2*(L2-M^2/L1)); %capacitor value

%frequency response from Laplace transform
w=[100e6:1e6:150e6]*2*pi;
V=(M/L1)*w0^2./(w0^2-w.^2);

%load impedance
ZL=(j*w0*(L1-M)+(1/(1/(j*w0*C)+j*w0*(L2-M))+1/(j*w0*M))^(-1));% at resonance
w=[100e6:1e5:150e6]*2*pi;
ZLw=(j*w.*(L1-M)+(1./(1./(j*w.*C)+j*w.*(L2-M))+1./(j*w.*M)).^(-1));
wp=1/sqrt(C*L2); %pole of impedance here
fp=wp/(2*pi)

figure;plot(w./(2*pi)*1e-6,abs(ZLw),'Linewidth',2);ylabel('|Z_L| (\Omega)');%'abs(input impedance)');
xlabel('Frequency (MHz)');axis([100 150 0 200]);%hold on;
title('Input Impedance');

