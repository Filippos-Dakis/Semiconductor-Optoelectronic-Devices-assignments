% Filippos Tzimkas-Dakis  UoC December 2020
% Final Project ------- Optoelectronics and Lasers
% Exercise 3 Quantum Well Gain
% Change the Eg_ =1.42 eV for GaAs and/or the QW width L
% also change the refractive index n = 3.52 (GaAS)
% variables with  X_ have units of eV
%%
tic
close all
clear all
clc

CD = cd;
addpath([CD,'\Functions_for_Consistency_Check'])


m0   = 9.11*10^(-31);     % electron mass (kg)
me   = 0.067*m0;          % e-mass in conduction band (kg)
mhh  = 0.45*m0;           % hh-mass (kg)
mlh  = 0.082*m0;          % lh-mass (kg)
mr_h = me*mhh/(me + mhh); % ehh reduced mass (kg)
mr_l = me*mlh/(me + mlh); % elh reduced mass (kg)

ev2joul = 1.60218e-19;    % converts eV to Joule
joul2ev = 6.242e+18;      % converts Joule to eV

kT_   = 0.026;               % (eV) Room temperature T=300K
kTT   = kT_*ev2joul;         % (joule) Room temperature T=300K
q     = 1.6*10^(-19);        % electric charge (Coulomb)
h     = 6.62607015*10^(-34); % Plancks Constant (J*s)
h_bar = h/2/pi;              % h-bar (J*s)

L  = 14 *10^(-9);          % QWell width (m) <<<<<-----------------

e0 = 8.85418*10^(-12);    % vacuum permitivity (Farad/m)
c  = 3*10^8;              % speed of light (m/s)
n  = 3.52;                % ref. index GaAs

Eg_  = 1.42;              % Energy bandgap (ev) <<< -------------
Ev_  = 0;                 % Energy reference level
Eg   = Eg_*ev2joul;       % Energy gap (joule)
Ee   =   (h_bar^2 *pi^2)/(2*me*L^2)*[1 2^2];  % Carrier QW energies (joule)
Ee_  = Ee*joul2ev;                            % Carrier QW energies (ev)
Ee__ = Ee_ + Eg_;
Ehh  = - (h_bar^2 *pi^2)/(2*mhh*L^2)*[1 2^2]; % HeavyHole QW energies (joule)
Ehh_ = Ehh*joul2ev;                           % HeavyHole QW energies (eV)
Elh  = - (h_bar^2 *pi^2)/(2*mlh*L^2)*[1 2^2]; % LightHole QW energies (joule)
Elh_ = Elh*joul2ev;                           % HeavyHole QW energies (eV)

const = pi*q^2 *h_bar/(2*n*e0*c*m0); %constant of g_max
Mthh  = [1/2 0];            % Mt^2 / M^2  heavy hole [TE TM]
Mtlh  = [1/6 2/3];          % Mt^2 / M^2  light hole [TE TM]
Mmo_  = 29;                 % (eV) for GaAs 2*(M^2)/m0
Mmo   = 29*ev2joul;         % (Joule) for GaAs 2*(M^2)/m0

E21_ehh = Eg + Ee - Ehh;    % eHH energy levels (joule)
E21_elh = Eg + Ee - Elh;    % eLH energy levels (joule)

r_ehh = mr_h/(pi*h_bar^2 *L)*[1 1]; % reduced density of states HH
r_elh = mr_l/(pi*h_bar^2 *L)*[1 1]; % reduced density of states LH

En21     = (0:0.001:3)*ev2joul;  % Energy vector (joule)

% Caltulating Max Gain for TE   ----------------------------------
% gain = 1/meters  for 1/cm  divide with 100 !!
Gmax_ehh_TE = const* 1 * r_ehh*Mmo *Mthh(1);  % HH TE polarization
Gmax_elh_TE = const* 1 * r_elh*Mmo *Mtlh(1);  % LH TE polarization
GMAX_TE     = [Gmax_ehh_TE(1) Gmax_elh_TE(1) Gmax_ehh_TE(2) Gmax_elh_TE(2)]; % STEPS
EE21_TE     = [E21_ehh(1) E21_elh(1) E21_ehh(2) E21_elh(2)];      % ROOTS
[~,ind1]    = min(abs(En21-E21_ehh(1))); % finding the spots of steps
[~,ind2]    = min(abs(En21-E21_elh(1))); % 
[~,ind3]    = min(abs(En21-E21_ehh(2))); % 
[~,ind4]    = min(abs(En21-E21_elh(2))); % 
% MaxGain envelope  1/m   TE  
Y_TE(1:ind1-1)          = 0;                          %
Y_TE(ind1:ind2-1)       = GMAX_TE(1);                 %
Y_TE(ind2:ind3-1)       = GMAX_TE(1) + GMAX_TE(2);    % Creating the Heaveside
Y_TE(ind3:ind4-1)       = GMAX_TE(1) + GMAX_TE(2) + GMAX_TE(3);% function
Y_TE(ind4:length(En21)) = sum(GMAX_TE);               %
Y_TE = (Y_TE/100)./En21;   % gain 1/cm
%------------------------------------------------------------
%------------------------------------------------------------
% Caltulating Max Gain for TM Polarization ------------------
% gain = 1/meters  for 1/cm  divide with 100 !!
Gmax_elh_TM = const* 1 * r_elh*Mmo *Mtlh(2);   % LH TE polarization  ONLY
GMAX_TM     = [Gmax_elh_TM(1) Gmax_elh_TM(2)]; % STEPS
EE21_TM     = [E21_elh(1) E21_elh(2)];         % ROOTS
[~,ind1]    = min(abs(En21-E21_elh(1))); % finding the spots of steps
[~,ind2]    = min(abs(En21-E21_elh(2))); % 
% MaxGain envelope  1/m   TM Polarization 
Y_TM(1:ind1-1)          = 0;                       %
Y_TM(ind1:ind2-1)       = GMAX_TM(1);              %
Y_TM(ind2:length(En21)) = GMAX_TM(1) + GMAX_TM(2); % Creating the Heaveside             
Y_TM = (Y_TM/100)./En21;   % gain 1/cm
% ----------------------------------------------------------


%% Carriers  and  Quasi-Fermi levels
E21_a   = 1.42 ;      % Minimun energy 
E21_b   = 3;          % Maximum energy
dE21    = 0.001 ;     % step
xLimits = [1 3];      % eV limits for plots

NP = [1:2:11 15 20 25 30];  % multiplyer of 10^12 carriers and holes
Fc = 0*NP;
Fv = 0*NP;

for i = 1:length(NP)
    % induced carriers and holes
    n_car  = NP(i)*10^12;      % 1/cm^2
    p_car  = n_car;            % 1/cm^2
    
    rhoe       = me/(pi*h_bar^2) ;    % Density of states 2D J^-1 m^-2
    rhoe_cm    = rhoe*ev2joul*10^-4 ; % 1/J 1/cm^2
    
    % density of states for the bound state
    he1        = @(x)rhoe_cm*(0-(-log(exp((Ee__(1) -x)/kT_)+1)*kT_ + Ee__(1)-x));
    % density of states for the second state
    he2        = @(x)rhoe_cm*(0-(-log(exp((Ee__(2)-x)/kT_)+1)*kT_  + Ee__(2)-x));
    
    h      = @(x)n_car - (he1(x)+he2(x)); % creates the function to be solved
    Efn   = fzero(h,1); % calculating Quasi-Fermi for CONDUCTION band
    Fc(i) = Efn; % stores the QF of Conduction band
    
    rhohh       = mhh/(pi*h_bar^2);     % Density of states 2D J^-1 m^-2
    rhohh_eV_cm = rhohh*ev2joul*10^-4 ; % 1/J 1/cm^2
    % density of states for the bound state
    jhh1        = @(x)rhohh_eV_cm*((log(exp((x-Ehh_(1))/kT_)+1)*kT_ + Ehh_(1) - x )-0);
    % density of states for the second state
    jhh2        = @(x)rhohh_eV_cm*((log(exp((x-Ehh_(2))/kT_)+1)*kT_ + Ehh_(2) - x )-0);
    
    rholh       = mlh/(pi*h_bar^2);       % Density of states 2D J^-1 m^-2                                                                                              
    rholh_eV_cm = rholh*ev2joul*10^-4;    % 1/J 1/cm^2
    % density of states for the bound state
    jlh1        = @(x)rholh_eV_cm*((log(exp((x-Elh_(1))/kT_)+1)*kT_ + Elh_(1) -x)-0);
    % density of states for the second state
    jlh2        = @(x)rholh_eV_cm*((log(exp((x-Elh_(2))/kT_)+1)*kT_ + Elh_(2)-x)-0);
    
    % creates the function to be solved
    j = @(x)p_car-(jhh1(x)+jhh2(x)+jlh1(x)+jlh2(x));
    
    Efp     = fzero(j,-1); % calculating Quasi-Fermi for Valance band
    Fv(i)   = Efp;  % stores the QF of valence band
    rhorelh = (mr_l/(pi*h_bar^2))/L; % reduced density of states for electron and HH J^-1 m^-2
    rhorehh = (mr_h/(pi*h_bar^2))/L; % reduced density of states for electron and LH J^-1 m^-2
    
    Mlh_TE = 29*ev2joul*m0/(2*6); % se eV->joule 1.6*10^-19 , 1/6 gia TE & LH
    Mhh_TE = 29*ev2joul*m0/(2*2);
    Mlh_TM = (2/3)*29*ev2joul*m0/2 ;% eV->joule, 2/3 for TM & LH
    Mhh_TM = 0;
    
    % initializing the functions of gain ---------------------------------
    gain_e_lh_TE = zeros((E21_b-E21_a)/dE21+1,1);
    gain_e_hh_TE = zeros((E21_b-E21_a)/dE21+1,1);
    gain_TE(:,i) = zeros((E21_b-E21_a)/dE21+1,1);
    
    gain_e_lh_TM = zeros((E21_b-E21_a)/dE21+1,1);
    gain_e_hh_TM = zeros((E21_b-E21_a)/dE21+1,1);
    gain_TM(:,i) = zeros((E21_b-E21_a)/dE21+1,1);
    % --------------------------------------------------------------------
    energy       = zeros((E21_b-E21_a)/dE21+1,1);
    
    jj = 0;
    for E21 = E21_a:dE21:E21_b
        jj = jj+1;
        
        % calculating the fermi distributions ----------
        E1lh  = Ev_ - (E21-Eg_)*mr_l/mlh;
        f1elh = 1/(exp((E1lh-Efp)/kT_)+1);
        
        E2lh  = Eg_+(E21-Eg_)*mr_l/me;
        f2elh = 1/(exp((E2lh-Efn)/kT_)+1);
        %-----------------------------------------------
        % checking in which are we ---------------------
        if E21 < Ee__(1)-Elh_(1)
            rhorelh_E21 = 0;
        elseif ( (E21 >= Ee__(1)-Elh_(1) ) && (E21 < Ee__(2)-Elh_(2)) )
            rhorelh_E21 = rhorelh;
        else
            rhorelh_E21 = 2*rhorelh; % when we go to the second state the
            % density of states doubles its value
        end
        %----------------------------------------------
        % Caculating the gain in     1/cm    -------------------
        % TE
        gain_e_lh_TE(jj) = (pi^1*q^2*h_bar*Mlh_TE*rhorelh_E21*(f2elh-f1elh)/(n*e0*c*m0^2*(E21*ev2joul)))/100;
        % TM
        gain_e_lh_TM(jj) = (pi^1*q^2*h_bar*Mlh_TM*rhorelh_E21*(f2elh-f1elh)/(n*e0*c*m0^2*(E21*ev2joul)))/100;
        %--------------------------------------------------------
        % calculating the fermi distributions ----------
        E1hh  = Ev_-(E21-Eg_)*mr_h/mhh;
        f1ehh = 1/(exp((E1hh-Efp)/kT_)+1);
        
        E2hh  = Eg_ +(E21-Eg_)*mr_h/me;
        f2ehh = 1/(exp((E2hh-Efn)/kT_)+1);
        %------------------------------------------------
        % checking in which are we ---------------------
        if E21 < Ee__(1)-Ehh_(1)
            rhorehh_E21=0;
        elseif ( (E21 >= Ee__(1)-Ehh_(1)) && (E21 < Ee__(2)-Ehh_(2)) )
            rhorehh_E21=rhorehh;
        else
            rhorehh_E21 = 2*rhorehh; % when we go to the second state the
            % density of states doubles its value
        end
        %------------------------------------------------
        % Caculating the gain in     1/cm    -------------------
        % TE
        gain_e_hh_TE(jj) = (pi^1*q^2*h_bar*Mhh_TE*rhorehh_E21*(f2ehh-f1ehh)/(n*e0*c*m0^2*(E21*ev2joul)))/100;
        % TM
        gain_e_hh_TM(jj) = (pi^1*q^2*h_bar*Mhh_TM*rhorehh_E21*(f2ehh-f1ehh)/(n*e0*c*m0^2*(E21*ev2joul)))/100;
        %--------------------------------------------------------
        % Calculating the TOTAL GAIN  TE and TM  in  1/cm
        gain_TE(jj,i) = gain_e_lh_TE(jj) + gain_e_hh_TE(jj);
        gain_TM(jj,i) = gain_e_lh_TM(jj) + gain_e_hh_TM(jj);
        
        energy(jj) = E21;  % storing for figures
    end
    
end

%% Prints !
fprintf('\n ========== PRINTS ARE HERE ===========')
fprintf('\n --------------- ROOTS ----------------\n')
fprintf(' Electron - HEAVY Hole (1) = %f eV \n',E21_ehh(1)*joul2ev)
fprintf('                       (2) = %f eV \n',E21_ehh(2)*joul2ev)
fprintf('\n Electron - LIGHT Hole (1) = %f eV \n',E21_elh(1)*joul2ev)
fprintf('                       (2) = %f eV \n',E21_elh(2)*joul2ev)
fprintf('\n ---------------------------------------\n')
fprintf('\n ------------ Quasi Fermi --------------\n')
for i = 1:length(NP)
    fprintf('\n For N =%d*10^(12), -> FC = %.3f eV \n',NP(i),Fc(i))
    fprintf('                    -> FV = %.3f eV \n', Fv(i))
end
fprintf('\n ======================================\n')

toc 
%% PLOTS !
% TE Plots
figure
% subplot(2,2,1)
subplot(1,2,1)
for i = 1:length(NP)
    plot(energy,gain_TE(:,i),'LineWidth',1.1);
    hold on
end
plot(En21*joul2ev,Y_TE,'LineWidth',1.5,'Color','r','LineStyle','-.')
plot(En21*joul2ev,-Y_TE,'LineWidth',1.5,'Color','r','LineStyle','-.')
plot(En21*joul2ev,0*Y_TE,'Color','k','LineStyle','-')
xlim([0.5 3]);
title('TE gain spectra')
xlabel('E_{21} ( eV )')
ylabel('Optical Gain (cm^{-1})')
xlim(xLimits)
hold off

% TM Plots
% figure
% subplot(2,2,2)
subplot(1,2,2)
for i = 1:length(NP)
    plot(energy,gain_TM(:,i),'LineWidth',1.1);
    hold on
end
plot(En21*joul2ev,Y_TM,'LineWidth',1.5,'Color','r','LineStyle','-.')
plot(En21*joul2ev,-Y_TM,'LineWidth',1.5,'Color','r','LineStyle','-.')
plot(En21*joul2ev,0*Y_TM,'Color','k','LineStyle','-')
xlim([0.5 3]);
title('TM gain spectra')
xlabel('E_{21} ( eV )')
ylabel('Optical Gain (cm^{-1})')
xlim(xLimits)
hold off 






