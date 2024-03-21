% Filippos Tzimkas-Dakis  UoC December 2020
% Final Project ------- Optoelectronics and Lasers
% Exercise 1 Tunable DBR Laser
%%
close all
clear all
clc 
tic

CD = cd;
addpath([CD,'\Functions_for_Consistency_Check'])


lambda0 = 1570*10^-9;  % center WL (meters)
n       = 3.4;         % refractive index
ng      = 3.8;         % group refractive index (Coldrens Table 2.1 pg 59)
hi      = 0.7;         % quantum efficiency
tau     = 3*10^(-9);   % carrier life time (sec)
WG_cs   = 0.2*3*10^(-12); % waveguide cross section (m2);
L_a     = 200*10^(-6); % active medium length (m)
L_p     = 200*10^(-6); % passive medium length (m)
L_g     = 200*10^(-6); % Bragg grating length (m)
k       = 100/(10^-2); % reflectivity/(m)
m2      = (ng -1)/(ng+1); % cleaved facet reflectvity
n_H     = 3.40;           % reflective index of InGaAsP @ 1.550
n_L     = 3.17;           % reflective index of InP @ 1.550
q       = 1.602*10^(-19); % electron charge Coulomb
dndN    = 10^(-21);       % cm^3
dndN    = dndN *10^(-6);  % m^3

% r     = (n_H - n_L)/(n_H + n_L); % reflectivity between two layers
LAMDA = lambda0/(2*n);         % n*LAMDA = lambda0 !
m     = L_g/LAMDA;             % number of periods 
Leff  = (1/2/k)*tanh(k*L_g);   % effective WL Coldrens (3.63)
% Leff_ = LAMDA/4/r;           % Approximation for high reflectivity  
dl0   = (lambda0^2)/(2*ng*(L_a + L_p + Leff)); % mode spacing (m)
Zg    = dndN * hi*tau/(q*WG_cs*L_g);  % everything from (3.77) except from I
                                      % 1/Ampere
Zp    = dndN * hi*tau/(q*WG_cs*L_p);  % same here 

Ig       = (0:0.0001:0.06);   % Bragg current (Ampere)
Delta_lg = lambda0*Zg*Ig/n;   % Bragg WL shift (nm)
%--------------------------------------------------------------------------
% mode shifting due to el. current in Bragg area
Delta_lm = (Delta_lg*n)*Leff/(ng *(L_a + L_p + Leff)); 
% New WL of lasing (only bragg current)
Delta_lasing1 = Delta_lm + round((Delta_lg - Delta_lm)/dl0)*dl0;
% -------------------------------------------------------------------------
%--------------------------------------------------------------------------
%%
Ip            = Ig*2.2636;        % Phase current (Ampere)
Delta_lp      = lambda0*Zp*Ip/n; % Phase WL shift (nm)
% mode shifting due to el. current in Bragg and Passive areas
Delta_lm2     = (Delta_lp*n*L_p +Delta_lg*n*Leff)/(ng *(L_a + L_p + Leff));
% New WL of lasing (Bragg and Phase current)
Delta_lasing2 = Delta_lm2 + round((Delta_lg - Delta_lm2)/dl0)*dl0;
%--------------------------------------------------------------------------
toc
%% Plotting results
figure
subplot(2,2,1)
plot(Ig*10^3,(lambda0 + Delta_lasing1)*10^9,'LineWidth',1.5,'Color','blue')
xlabel('Current I_g ( mA )')
ylabel('Lasing Wavelength ( nm )')
title('Only I_g')
hold off

subplot(2,2,2)
plot(Ig*10^3,Delta_lg*10^9,'LineWidth',1.5,'Color','blue')
hold on
plot(Ig*10^3,Delta_lm*10^9,'LineWidth',1.5,...
    'LineStyle','-.','Color','r')
xlabel('Current I_g ( mA )')
ylabel('\Delta _\lambda ( nm )')
legend('\Delta\lambda _g','\Delta\lambda_m')
title('Only I_g')
hold off

subplot(2,2,3)
plot(Ig*10^3,(lambda0 + Delta_lasing2)*10^9,'LineWidth',1.5,'Color','blue')
xlabel('Current I_g ( mA )')
ylabel('Lasing Wavelength ( nm )')
title('I_g and I_p')
hold off

subplot(2,2,4)
plot(Ig*10^3,Delta_lg*10^9,'LineWidth',1.5,'Color','blue')
ylabel('\Delta _\lambda ( nm )')
xlabel('Current I_g ( mA )')
hold on
plot(Ig*10^3,Delta_lm2*10^9,'LineWidth',1.5,...
    'LineStyle','-.','Color','r')
legend('\Delta\lambda _g','\Delta\lambda_m')
title('I_g and I_p')
% ax1 = gca; % current axes
% ax1_pos = ax1.Position; % position of first axes
% ax2 = axes('Position',ax1_pos,...
%     'XAxisLocation','top',...
%     'YAxisLocation','right',...
%     'Color','none');
% line(Ip*10^3,Delta_lm2*10^9,'LineWidth',1.5,...
%     'Parent',ax2,'LineStyle','-','Color','r')
% xlim([0 max(Ip)*10^3])6.592 "
%xlabel('Current I_g ( mA )')
%%
figure
plot(Ig*10^3,(lambda0 + Delta_lasing1)*10^9,'LineWidth',1.5,'Color','blue')
hold on
plot(Ig*10^3,(lambda0 + Delta_lasing2)*10^9,'LineWidth',1.5,'Color','red')
xlabel('Current I_g ( mA )')
ylabel('Lasing Wavelength ( nm )')
hold off
%% Plotting both cases (a) and (b) in common graph
[~,ind] = min(abs(Ip-0.05));
figure
% subplot(2,2,1)
yyaxis left
plot(Ig(1:ind)*10^3,(lambda0 + Delta_lasing1(1:ind))*10^9,...
                        'LineWidth',1.5,'Color','blue')
hold on
plot(Ig(1:ind)*10^3,(lambda0 + Delta_lasing2(1:ind))*10^9,...
                       'LineWidth',1.5,'Color','red')
xlabel('Current I_g ( mA )')
ylabel('Lasing Wavelength ( nm )')
ylim((lambda0+[0,Delta_lasing2(ind)*1.01])*10^9)
xlim([0 Ig(ind)]*10^3)
yyaxis right
plot(Ig(1:ind)*10^3, Ip(1:ind)*10^3,'LineWidth',1.5,'Color','k')
ylabel('Phase Current ( mA )')
ylim([0 60])
xlabel('Bragg Current I_p ( mA )')
plt = gca;
plt.YAxis(1).Color = 'k'; 
plt.YAxis(2).Color = 'k'; 
                      


                      
                      