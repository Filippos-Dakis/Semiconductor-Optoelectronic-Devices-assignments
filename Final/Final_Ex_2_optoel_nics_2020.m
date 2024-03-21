% Filippos Tzimkas-Dakis  UoC December 2020
% Final Homework ------- Optoelectronics and Lasers
% Exercise 2 Bragg reflector
%%
close all
clear all
clc 
tic

CD = cd;
addpath([CD,'\Functions_for_Consistency_Check'])

n_h = 2.1; % High Refractive index
n_l = 1.7; % Low Refractive index
n_a = 1.5;  % Refr. Index from the left
n_b = 1.9;  % Refr. Index from the right
%%
Threshold = 0.99; % reflectivity at center WL
Gamma     = 0;    % Reflection coefficient Gamma (transmission lines, microwaves etc)
                  % it is connected with r, the reflection coefficient of the EM fiel
R         = Gamma^2; % Power Reflection Coefficient R^2
N         = 1;    % intialize N
% Calculating the periodes needed for a specified reflectivity
while (R<= Threshold )
    Gamma = (n_a*(n_h/n_l)^(2*N)-n_b)/(n_a*(n_h/n_l)^(2*N)+n_b);
    R     = Gamma^2;
    N = N +1;
end
N = N-1; % Total number of periods needed
fprintf('\n The periods needed are %d \n',N)
fprintf(' The Reflectivity is R = %f \n \n',R)
%%
c       = 3*10^8;               % speed of light (m/s)
lambda0 = 600*10^-9;            % center wavelength (meters)
k0      = 2*pi/lambda0;         % wave-number at center wavelength 
lambda  = (300:0.1:800)*10^-9; % Wavelength span for sweep
d_lambda = lambda(2) - lambda(1);
freq    = linspace(0,3/lambda0,2000);% frequebcy span for sweep
df      = freq(2) - freq(1);
k       = 2*pi./lambda;         % Wavevector in free space
d_h     = (1/n_h)*(lambda0/4);  % physical Q-length for High-index layer
d_l     = (1/n_l)*(lambda0/4);  % physical Q-length for Low-index layer
k_h     = k*n_h;                % Wavevector in High index layer
k_l     = k*n_l;                % Wavevector in Low indez layer
heta    = 377;                  % free space characteristic impedance
Z_h     = heta/n_h;             % High-index characteristic impedance
Z_l     = heta/n_l;             % Low-index characteristic impedance
M       = 2*N;
Z       = 0*ones(M+1,length(lambda));

Zf(M+1,1:length(freq)) = (heta/n_b)*ones(1,length(freq));    % Char. Impedance at the Right End (material b) FREQ
Z(M+1,:)  = heta/n_b;    % Char. Impedance at the Right End (material b) WVLG
Z_a      = heta/n_a;    % Char. Impedance at the Left End (material a)
% Setting the parameters for the periodic stack
for i = M:-2:2
    ZZ(i)    = Z_h;
    d(i)     = d_h;
    kk(i,:)  = k_h; % k-number in High index layer for wavelengthe sweep
    kkf(i,:) = 2*pi*freq*n_h; % k-number in High index layer for frequency sweep
    
    ZZ(i-1)    = Z_l;
    d(i-1)     = d_l;
    kk(i-1,:)  = k_l;% k-number in Low index layer for wavelengthe sweep
    kkf(i-1,:) = 2*pi*freq*n_l; % k-number in Low index layer for frequency sweep
end
%%
for q = M+1:-1:2
     % calculating characteristic impedance for wavelength sweep
     Z(q-1,:) = ZZ(q-1)*(Z(q,:) + 1i*ZZ(q-1)*tan(kk(q-1,:)*d(q-1)))./(ZZ(q-1) + 1i*Z(q,:).*tan(kk(q-1,:)*d(q-1)));
     % calculating characteristic impedance for frequency sweep
     Zf(q-1,:)=ZZ(q-1)*(Zf(q,:) + 1i*ZZ(q-1)*tan(kkf(q-1,:)*d(q-1)))./(ZZ(q-1) + 1i*Zf(q,:).*tan(kkf(q-1,:)*d(q-1)));
end
%%
% Reflection coefficient Gamma (transmission lines, microwaves etc)
% it is connected with r, the reflection coefficient of the EM field
G  = (Z(1,:)- Z_a)./(Z(1,:) + Z_a);
% Reflection coefficient Gamma for frequency sweep
Gf = (Zf(1,:)- Z_a)./(Zf(1,:) + Z_a);
% Power Reflection Coefficient R^2
RR = abs(G).^2;
% Power Transmission Coefficient T^2
TT = 1 - RR;
% Power Reflection Coefficient R^2, frequency sweep
RRf = abs(Gf).^2;
% Power Transmission Coefficient T^2, frequency sweep
TTf = 1-RRf;
RRf_dB = 10*log10(RRf);
RR_dB = 10*log10(RR);
% Calculating the 3-dB Bandwidth of the Reflector/Mirror
BW_f  = fwhm_dakis(RRf')*df*c;  % BW in frequency
BW_wl = fwhm_dakis(RR')*d_lambda*c/(lambda0^2); % BW in wavelength

toc
%% Plotting  Results 
figure
subplot(2,2,1)
plot(lambda*10^9,RR,'LineWidth',1.2)
hold on
xlim([min(lambda) max(lambda)]*10^9);
set(gca,'FontSize',14);
ylabel(' R ')
xlabel(' Wavelength  ( nm ) ')
hold off

subplot(2,2,3)
plot(lambda*10^9,RR_dB,'LineWidth',1.2,'Color','blue')
xlim([min(lambda) max(lambda)]*10^9);
set(gca,'FontSize',14);
ylim([-30,0.5])
ylabel(' R  ( dB ) ')
xlabel(' Wavelength  ( nm ) ')
hold off

subplot(2,2,2)
plot(freq*lambda0,RRf,'LineWidth',1.2)
xlim([min(freq*lambda0) max(freq*lambda0)]);
set(gca,'FontSize',14);
ylabel(' R ')
xlabel(' f/f_0 ')
hold off

subplot(2,2,4)
plot(freq*lambda0,RRf_dB,'LineWidth',1.2,'Color','b')
xlim([min(freq*lambda0) max(freq*lambda0)]);
ylim([-30 0.5])
set(gca,'FontSize',14);
ylabel(' R  ( dB )')
xlabel(' f/f_0 ')
hold off
%% Online code cross-Check 
na = n_a; 
nb = n_b; 
nH = n_h; 
nL = n_l; % refractive indices

LH = d_h*n_h/lambda0;
LL = d_l*n_l/lambda0; % optical thicknesses in units of ?0
la0 = 600; % ?0 in units of nm
rho = (nH-nL)/(nH+nL); % reflection coefficient ?
la2 = pi*(LL+LH)*1/acos(rho) * la0; % right bandedge
la1 = pi*(LL+LH)*1/acos(-rho) * la0; % left bandedge
Dla = la2-la1; % bandwidth
N = M/2; % number of bilayers
n = [na, repmat([nL,nH], 1, N), nb]; % indices for the layers A|H(LH)N|G
L = [repmat([LL,LH], 1, N)]; % lengths of the layers H(LH)N
la = linspace(100,1000,2001); % plotting range is 300 ? ? ? 800 nm
Gla = abs(multidiel(n,L,la/la0)).^2; % reflectance as a function of ?
figure; plot(la,Gla);



