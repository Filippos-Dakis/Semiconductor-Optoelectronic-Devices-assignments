% Filippos Tzimkas-Dakis  UoC October 2020
% Homework 1 ------- Optoelectronics and Lasers


% Solving the square welll potential
%------ Even solutions -------------
eV2J      = 1.60218e-19;            % converts eVolts to Joules
J2eV      = 6.241509e+18;           % converts Joules to eVolts
% c       = 3*1e8;                  % speed of light (m/s)
Eg_GA     = 1.51;                   % E-Bandgap @ GaAs(eV)
Eg_AGA    = 1.91;                   % E-Bandgap @ AlGaAs(eV)
Delta_Eg  = Eg_AGA -Eg_GA;
Eg        = Eg_GA;
me        = 9.1093837015*1e-31;          % electron mass in (Kg)
h_eV      = 6.582119569*1e-16;           % h-bar in (eV*s)
h_J       = 1.054571817*1e-34;           % h-bar in (J*s)
L         = (1:0.02:50)*1e-9;            % creates the a vector with widths
% h_ev^2 /me
%% Pure Mathematical way for ee1
m    = 0.07*me;               % effective mass of the electron
Vo   = (2/3)*Delta_Eg;        % well's potential in (eV)
Vo   = Vo*eV2J;
Uo   = 2*m*Vo/h_J^2;          % defines a new normalized potential value
y    = sqrt(Uo)*L/2;

Eee     = 0.*L;
options = optimset('Display','off'); %meaningless input for fsolve
out     = 0.1;

for i = 1:length(L)
    % using "fsolve()" to solve the trascedental equation
    out  = fsolve(@(t)tan(t) - sqrt(y(i)^2 - t^2)/t,out,options);
    % Eee(i) stores the g.s. energy for every L
    Eee(i) = (2*out/L(i))^2*(h_J^2 /2/ m)*J2eV;
end
%% Pure Mathematical way for hh1
m    = 0.4*me;               % effective mass of the electron
Vo   = (1/3)*Delta_Eg;       % well's potential in (eV)
Vo   = Vo*eV2J;
Uo   = 2*m*Vo/h_J^2;         % defines a new normalized potential value
y    = sqrt(Uo)*L/2;

Ehh     = 0.*L;
options = optimset('Display','off'); % meaningless input for fsolve
out     = 0.1;

for i = 1:length(L)
    % fsolve, solves the trascedental equation
    out  = fsolve(@(t)tan(t) - sqrt(y(i)^2 - t^2)/t,out,options); %solver
    % Ehh(i) stores the g.s. energy for every L
    Ehh(i) = (2*out/L(i))^2*(h_J^2 /2/ m)*J2eV;
end

%%
close all
E_tot = Eg + Ehh + Eee;

figure(1)                            % Ploting results
plot(L*10^9,Eee,'LineWidth',1.5)
hold on
plot(L*10^9,Ehh,'LineWidth',1.5)
plot(L*10^9,E_tot,'LineWidth',1.5)
axis on
xlabel('Well`s width (nm)')
ylabel('Energy (eV)')
legend('ee_1','hh_1','E_g + ee_1 + hh_1')

E1 = 1.55174; % eV @ 799 nm
E2 = 1.66422; % ev @ 745 nm
[~,j1] = min(abs(Eg + Ehh+Eee - E1));
L(j1)                                % returns the width L_1
[~,j2] = min(abs(Eg + Ehh+Eee -E2));
L(j2)                                % returns the width L_2
% 1.55174 ev  @ 799
% 1.55564 ev  @ 797
% 1.66422 ev  @ 745


