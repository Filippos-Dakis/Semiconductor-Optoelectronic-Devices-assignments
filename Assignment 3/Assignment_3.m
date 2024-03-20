% Filippos Tzimkas-Dakis  UoC November 2020
% Homework 3 ------- Optoelectronics and Lasers
% Script for the Bonus exercise of homework 3
%%
clc
clear all
close all

c     = 3*10^10;      % speed of light in vacuum (cm/s)
lamda = 1.3*10^(-4); % um = 10^-6 m = 10^-4 cm
G     = 0.2;          % confinement (clear)
b_sp  = 10^-5;        % clear
tp    = 1*10^(-12);   % photon lifetime (s)
B     = 10^(-10);     % cm^3/s
H     = 200*10^(-7);  % active region thickness (cm)
L     = 500*10^(-4);  % cavity length (cm)
W     = 50*10^(-4);   % width of metal contacts (cm)
Vp    = H*L*W;        % Cavity Volume (cm^-3)
V     = G*Vp;         
hi    = 1;            % injection efficiency (clear)
a     = 2.5*10^(-16); % cm^2
N_tr  = 2*10^(18);    % 1/(cm^3);
n1    = 3.2;          % refractive index (clear)
dn_dl = -1*10^(+4);   % dn/dlamda = -1/(um)
n     = n1 - lamda*dn_dl; % refractive index for InP at 1.55 um
ug    = c/n;          % group velocity
q     = 1.602*10^-19; % electric charge
T     = 300;          % Temperature (Kelvin)
I     = (0:1:1000)*10^-3; % Ampere
% intialize variables 
N1    = NaN*ones(1,length(I));
N2    = N1;
N3    = N2;

% define new parameters to help the solver
K  = G*b_sp*tp*a*ug; 
P  = G*ug*tp*a;
Z  = hi*I/q/V;
%%
syms X
for i = 1:length(I)
    % characteristic equation
    EQN = B*(K-P)*X^3 + B*(1 + P*N_tr - K*N_tr)*X^2 + Z(i)*P*X - (1+P*N_tr)*Z(i)==0;
    % we use VPAsolve because fsolse couldn't help
    temp = vpasolve(EQN,X);
    N1(i) = temp(1);    % solution separation 
    N2(i) = temp(2);    % we are interested in this
    N3(i) = temp(3);
end

[~,index] = min(abs(N2-N3)); % finds threshold index
Ith      = I(index);         % threshold current 
Nth      = N2(index);        % thhreshold carriers
N_nrmlzd = N2/Nth;           % normalized carriers density
%%
% photon density
Np = (tp*b_sp*G*B*N2.^2)./(1-tp*G*ug*a*(N2-N_tr));

f1 = figure; % plot for photon density
semilogy(I*10^3,Np,'LineWidth',2)
xlabel('Device Current (mA)')
ylabel('N_p (cm^{-3}) ')
title('Photon Density')
hold off
f1.Units = 'normalized';
f1.OuterPosition = [0.5398    0.5200    0.3367    0.4756];

f2 = figure; % plot for carrier density
plot(I*10^3,N_nrmlzd,'LineWidth',2)
xlabel('Device Current (mA)')
ylabel('N (cm^{-3}) ')
legend('Normalized')
ylim([-0 1.1])
title('Carrier Density')
hold off
f2.Units = 'normalized';
f2.OuterPosition = [0.0988    0.5188    0.3367    0.4756];

%%
f3 = figure; % all solutions together
plot(I*10^3,N1,'LineWidth',2,'Color','b')
hold on
plot(I*10^3,N2,'LineWidth',2,'Color','black')
plot(I*10^3,N3,'LineWidth',2,'Color','red')
xlabel('Device Current (mA)')
ylabel('N (cm^{-3}) ')
legend('N_1','N_2','N_3')
hold off
f3.Units = 'normalized';
f3.OuterPosition = [0.3293    0.0519    0.3367    0.4756];
%%
options = optimset('Display','off');
    fun1 = @(x)B*(K-P)*(Nth^3)*x^3 + B*(1 + P*N_tr - K*N_tr)*(Nth^2)*x^2 + Z(1)*P*Nth*x - (1+P*N_tr)*Z(1);
    % we use VPAsolve because fsolse couldn't help
    temp = fsolve(fun1,eps,options);
for i = 2:100
    % characteristic equation
    fun1 = @(x)B*(K-P)*(Nth^3)*x^3 + B*(1 + P*N_tr - K*N_tr)*(Nth^2)*x^2 + Z(i)*P*Nth*x - (1+P*N_tr)*Z(i);
    % we use VPAsolve because fsolse couldn't help
    temp = fsolve(fun1,-eps,options);
    N1(i) = temp;    % solution separation 

end





