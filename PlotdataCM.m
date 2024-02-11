close all;clear;
% If new data is used for this file, may need to change code based on
% \DeltaTau or organization of the data
%% Load in Transfer Function Data
% Plane Wave Simulation: data given on 2023_07_02
load('Efsweep2'); %magnitude
load('Efphase2'); %phase
Efsweep(:,2)=Efsweep2;Efphase(:,2)=Efphase2; 
Efsweep(:,1)=0:20:200;Efphase(:,1)=0:20:200;

%Reciprocity lab experiment:
load('RecExpReal1'); %real
load('RecExpImag1'); %imaginary
load('RecExpReal2'); 
load('RecExpImag2'); 
load('RecExpReal3'); 
load('RecExpImag3'); 
load('RecExpRealBG');
load('RecExpImagBG');
BGR = RecExpRealBG + 1j*RecExpImagBG;
if 0 %if subtraction background for the reciprocity experiments
    RecExpM1=(abs(RecExpReal1+1j*RecExpImag1 - BGR));
    RecExpP1=(angle(RecExpReal1+1j*RecExpImag1 - BGR));
    RecExpM2=(abs(RecExpReal2+1j*RecExpImag2 - BGR));
    RecExpP2=(angle(RecExpReal2+1j*RecExpImag2 - BGR));
    RecExpM3=(abs(RecExpReal3+1j*RecExpImag3 - BGR));
    RecExpP3=(angle(RecExpReal3+1j*RecExpImag3 - BGR));
else
    RecExpM1=(abs(RecExpReal1+1j*RecExpImag1));
    RecExpP1=(angle(RecExpReal1+1j*RecExpImag1));
    RecExpM2=(abs(RecExpReal2+1j*RecExpImag2));
    RecExpP2=(angle(RecExpReal2+1j*RecExpImag2));
    RecExpM3=(abs(RecExpReal3+1j*RecExpImag3));
    RecExpP3=(angle(RecExpReal3+1j*RecExpImag3));
end


%Reciprocity simulation
load('rTFSimple.mat'); %rTF - tansfer function for simple lead and voltage source near the tip
load('rTFSimplePhase.mat'); 
rTFphase=rTFphase*pi/180;

%Loop coil excitor simulation - Background already subtracted
EfPix(:,1)=0:5:200;EfPixPhase(:,1)=0:5:200;
load('Eypix'); %no extrapolation - may want to extrapolate for the last data point tau =200mm
EfPix(:,2)=Eypix(1,:);EfPixPhase(:,2)=Eypix(2,:);

%Experimental data - 3 measurements and 1 background
tas=5:45;% have 49 measurements, extra 4*5=20mm each side
load('TF20cm_BG.mat'); 
BG=Value_of_Ch1_R_0s(tas)+1i*Value_of_Ch1_I_0s(tas);
load('TF20cm_meas1.mat');
m1=Value_of_Ch1_R_0s(tas)+1i*Value_of_Ch1_I_0s(tas)-BG;
load('TF20cm_meas2.mat');
m2=Value_of_Ch1_R_0s(tas)+1i*Value_of_Ch1_I_0s(tas)-BG;
load('TF20cm_meas3.mat');
m3=Value_of_Ch1_R_0s(tas)+1i*Value_of_Ch1_I_0s(tas)-BG;
tau=0:5:200; 


%% Normalize the transfer functions
%reciprocity experiment
RecExpM1=RecExpM1/(trapz(RecExpM1)*5); %normalize
RecExpM2=RecExpM2/(trapz(RecExpM2)*5);
RecExpM3=RecExpM3/(trapz(RecExpM3)*5);
allMR = [RecExpM1 RecExpM2 RecExpM3]; 
AvgMR=mean(allMR,2); % average (center) of magnitude data to plot
semMR = std(allMR,0,2)/sqrt(3); % SEM for magnitude
RecExpP1=unwrap(RecExpP1);RecExpP1=RecExpP1-RecExpP1(1);
RecExpP2=unwrap(RecExpP2);RecExpP2=RecExpP2-RecExpP2(1);
RecExpP3=unwrap(RecExpP3);RecExpP3=RecExpP3-RecExpP3(1);
allPR = [RecExpP1 RecExpP2 RecExpP3]; 
AvgPR=mean(allPR,2);
semPR = std(allPR,0,2)/sqrt(3);
RecExpTau=0:5:200;

%for piecewise lab experiment
Dtau=5;
mm1=abs(m1)/(trapz(abs(m1))*Dtau);
mm2=abs(m2)/(trapz(abs(m2))*Dtau);
mm3=abs(m3)/(trapz(abs(m3))*Dtau);
%Get average and standard error of the mean (SEM) 
allM = [mm1 mm2 mm3]; 
AvgM=mean(allM,2); % average (center) of magnitude data to plot
semM = std(allM,0,2)/sqrt(3); % SEM for magnitude
m1a=unwrap(angle(m1)); m2a=unwrap(angle(m2));m3a=unwrap(angle(m3));
m1a=m1a-m1a(1); 
m2a=m2a-m2a(1);
m3a=m3a-m3a(1);
allP = [m1a m2a m3a];
AvgP=mean(allP,2);
semP = std(allP,0,2)/sqrt(3);

%for loop coil simulation
EfPix(:,2)=EfPix(:,2)/(trapz(EfPix(:,2))*Dtau);
EfPixPhase(:,2)=unwrap(EfPixPhase(:,2));
EfPixPhase(:,2)=EfPixPhase(:,2)-EfPixPhase(1,2);

%for plane wave simulation
DtauPW=200/(length(Efsweep(:,2))-1);
Efsweep(:,2)=Efsweep(:,2)/(trapz(Efsweep(:,2))*DtauPW);
Efphase(:,2)=unwrap(Efphase(:,2));
Efphase(:,2)=Efphase(:,2)-Efphase(1,2);

%for reciprocity simulation
rTFphase=unwrap(rTFphase);rTFphase=rTFphase-rTFphase(1);
DeltaTau=200/(length(rTF)-1);%normalize for 20cm long lead
rTF=rTF/(trapz(rTF)*DeltaTau);
tauP=[0:DeltaTau:200];

%% Plot Transfer Functions
figure;
subplot(2,1,1);
plot(Efsweep(:,1),Efsweep(:,2),'b','Linewidth',2);
ylabel('Normalized Magnitude |S|');title('Transfer functions');
hold on;
plot(EfPix(:,1),EfPix(:,2),'--','color','r');hold on;
plot(tauP,rTF,'--','Linewidth',2,'color','m'); hold on;
errorbar(tau,AvgM,semM); %plot the nice plots for SEM 
%plot(RecExpTau,RecExpM,':','Linewidth',2);
errorbar(RecExpTau,AvgMR,semMR); 
legend('Plane Wave Simulation','Excitor Simulation','Reciprocity Simulation','Excitor Experiment','Reciprocity Experiment');

subplot(2,1,2);
plot(Efphase(:,1),Efphase(:,2),'b','Linewidth',2);ylabel('Phase (radians) \angleS');xlabel('\tau (mm)'); hold on;
plot(EfPixPhase(:,1),EfPixPhase(:,2),'--','color','r');hold on;
plot(tauP,rTFphase,'--','Linewidth',2,'color','m'); hold on;
errorbar(tau,AvgP,semP);hold on;
%plot(RecExpTau,RecExpP,':','Linewidth',2);
errorbar(RecExpTau,AvgPR,semPR); 

%% Calculate Percent Errors
A=Efsweep(:,2);Ap=Efphase(:,2); % plane wave - A
B=EfPix(:,2);Bp=EfPixPhase(:,2); % loop coil - B
C=AvgM;Cp=AvgP; % loop coil experiment       - C
D=rTF;Dp=rTFphase; % reciprocity simulations - D
% E=RecExpM.';Ep=RecExpP.';%reciprocity experiment - E
E=AvgMR.';Ep=AvgPR.';%reciprocity experiment - E


dtA=200/(length(A)-1);dtB=200/(length(B)-1);
dtC=200/(length(C)-1);dtD=200/(length(D)-1);
dtE=5; 

tA=[0:dtA:200];tB=[0:dtB:200];tC=[0:dtC:200];tD=[0:dtD:200];
t=[0:5:200];
%interpolate everything to [0:5:200]
A=interp1(tA,A,t);
B=interp1(tB,B,t);
C=interp1(tC,C,t);
D=interp1(tD,D,t);
%E=interp1(tE,E,t);

Ap=interp1(tA,Ap,t);
Bp=interp1(tB,Bp,t);
Cp=interp1(tC,Cp,t);
Dp=interp1(tD,Dp,t);
%Ep=interp1(tE,Ep,t);

% experiment vs plane wave sim, pix sim, and reci sim
AllData=[A; B; D]; %plane wave, loop coil, reciprocity
AllDatap=[Ap; Bp; Dp]; %plane wave, loop coil, reciprocity

for i=1:3
    ErrorC(i,:)=abs((AllData(i,:)-C)./C)*100;
    ErrorCp(i,:)=abs((AllDatap(i,:)-Cp)./Cp)*100;
    ErrorE(i,:)=abs((AllData(i,:)-E)./E)*100;
    ErrorEp(i,:)=abs((AllDatap(i,:)-Ep)./Ep)*100;
end
xmin=0; xmax=200; ymin=0; ymax=500;

figure;
subplot(2,1,1);
plot(tau,ErrorC(1,:),'b','Linewidth',2);hold on;ylabel('% Error in |S|');
plot(tau,ErrorC(2,:),'--','color','r');hold on;ylabel('% Error in |S|');
plot(tau,ErrorC(3,:),'--','color','m','Linewidth',2);hold on;ylabel('% Error in |S|');
title('% Error Compared to Piecewise Method');
axis([xmin xmax ymin ymax]);

subplot(2,1,2);
plot(tau,ErrorCp(1,:),'b','Linewidth',2);hold on;ylabel('% Error in \angleS');xlabel('\tau (mm)');
plot(tau,ErrorCp(2,:),'--','color','r');hold on;ylabel('% Error in \angleS');xlabel('\tau (mm)');
plot(tau,ErrorCp(3,:),'--','color','m','Linewidth',2);hold on;ylabel('% Error in \angleS');xlabel('\tau (mm)');
legend('Plane Wave Simulation','Excitor Simulation','Reciprocity Simulation');
axis([xmin xmax ymin ymax]);

figure;
subplot(2,1,1);
plot(tau,ErrorE(1,:),'b','Linewidth',2);hold on;ylabel('% Error in |S|');
plot(tau,ErrorE(2,:),'--','color','r');hold on;ylabel('% Error in |S|');
plot(tau,ErrorE(3,:),'--','color','m','Linewidth',2);hold on;ylabel('% Error in |S|');
title('% Error Compared to Reciprocity Method');
axis([xmin xmax ymin ymax]);

subplot(2,1,2);
plot(tau,ErrorEp(1,:),'b','Linewidth',2);hold on;ylabel('% Error in \angleS');xlabel('\tau (mm)');
plot(tau,ErrorEp(2,:),'--','color','r');hold on;ylabel('% Error in \angleS');xlabel('\tau (mm)');
plot(tau,ErrorEp(3,:),'--','color','m','Linewidth',2);hold on;ylabel('% Error in \angleS');xlabel('\tau (mm)');
legend('Plane Wave Simulation','Excitor Simulation','Reciprocity Simulation');
axis([xmin xmax ymin ymax]);

%get average percent errors for insulated regions
AvgC=mean(ErrorC(:,3:end-3),2)
AvgCp=mean(ErrorCp(:,3:end-3),2)
AvgE=mean(ErrorE(:,3:end-3),2)
AvgEp=mean(ErrorEp(:,3:end-3),2)

%get average percent errors for lead tips
AvgCt=mean([ErrorC(:,1:2)  ErrorC(:,end-2:end)],2)
AvgCpt=mean([ErrorCp(:,2)  ErrorCp(:,end-2:end)],2)
AvgEt=mean([ErrorE(:,1:2)  ErrorE(:,end-2:end)],2)
AvgEpt=mean([ErrorEp(:,2)  ErrorEp(:,end-2:end)],2)

%overall percent error
AvgCo=mean(ErrorC,2)
AvgCpo=mean(ErrorCp(:,2:end),2)
AvgEo=mean(ErrorE(:,1:end-1),2) 
AvgEpo=mean(ErrorEp(:,2:end-1),2)


%find correlations
AA=[A; B; C; D; E].';
AAA=corrcoef(AA);
figure;imagesc(AAA);xticks([]);yticks([]);colorbar;
%xticklabels({' ' 'Loop Coil Simulation'});


%% PiX excitor 128HP data
%reflection coefficient of simulated Excitor with 40.5pF capacitor in
%sim4life 
load('ExcitorAbsreflectioncoeff'); % Simulation reflection coefficient
load('ExperimentRefl'); % reflection coefficient from experiment (from datasheet in Zurick Med Tech box)

figure;
plot(Gamma(:,1)*1e3,Gamma(:,2),'--'); hold on;
ylabel('Reflection Coefficient (dB)');xlabel('Frequency (MHz)');
title('Reflection coefficient of excitor');axis([100 180 -18 -4]);
plot(ExperimentRefl(:,1),ExperimentRefl(:,2));
legend('FDTD simulation','Experiment with network analyzer');

load('CapV'); % The voltage at the capactior for the same Sim4life simulation as reflection coeff above
figure;
plot(CapV(:,1)*1e3,CapV(:,2)); title('Capacitor voltage'); ylabel('Voltage Magnitude');
xlabel('Frequency (MHz)');axis([70 220 0 1e-10]);
% find minimum for simulation reflection coefficient and maximum for
% capacitor voltage
vmax = max(CapV(:,2));
index1 = find(CapV(:,2) >= vmax, 1, 'first');

load('Exat10mmexcitor');%data that shows tangential electric field 10mm away from capacitor
figure;
plot(Exat10mm(:,1),Exat10mm(:,2)); hold on;
x1=[-70:.01:70];
y1=[zeros(1,6750) ones(1,length(x1)-6750*2) zeros(1,6750)];%[zeros(1,67.5*2) ones(1,11) zeros(1,67.5*2)];
plot(x1,y1,'--');
xlabel('Distance from part of lead with shortest distance to capacitor (mm)');
ylabel('RMS of E_t_a_n (V/m)');
legend('Loop Coil Excitor','Plane Wave Sources');title('Electric field tangent to lead');
%find FWHM of the PiX excitation
y=Exat10mm(:,2);x=Exat10mm(:,1);
halfMax = max(y) / 2;
% Find where the data first drops below half the max.
index1 = find(y >= halfMax, 1, 'first');
% Find where the data last rises above half the max.
index2 = find(y >= halfMax, 1, 'last');
fwhm = index2-index1 + 1; % FWHM in indexes
% if you have an x vector
fwhmx = x(index2) - x(index1);


