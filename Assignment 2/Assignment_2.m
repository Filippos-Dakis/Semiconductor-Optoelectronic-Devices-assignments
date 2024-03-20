% Filippos Tzimkas-Dakis  UoC October 2020
% Homework 2 ------- Optoelectronics and Lasers
% This is a script to calculate the waveguided eigenmodes TM of a field in
% a dielectric slab waveguide.
% Also, this script calculates all the transmision constants \beta for
% specific characteristics (refractive indexes n1,n2,n3, width/thickness
%h, frequency or wave length f,\lambda).
%% Solving for TE modes
close all
clear all
clc
c  = 3*10^8; % speed of light in vacuum
n1 = 2.6;    % refractive index of the waveguide layer / core
n2 = 2.4;      % refractive inde of the substrate
n3 = 1;     % refractive inde of the cladding layer

d      = 500*1e-9;     % thickness of the eaveguide layer (meters)
lamda0 = 350*1e-9;     % wavelength in vacuum
k0     = 2*pi/lamda0;  % wave-number in vacuum
k1     = n1 *k0;       % wave-number in core
k2     = n2 *k0;       % wave-number in substrate
k3     = n3 *k0;       % wave-number in cladding layer

xi = (n2^2-n3^2)/(n1^2-n2^2); % APDWG "asymmetry" factor (0==SPDWG)
V = k0*d/2*sqrt( n1^2-n2^2 ); % V-number, same definition as for SPDWG
VasosTE = atan( sqrt(xi) )/2; % V-number "asymmetry-offset" for TE-modes
TE_modes = ceil( (V-VasosTE)/(pi/2) ); % Number of modes expected

%========== The waveguide looks like this ==========
%|                                      |  y?
%|     n3 , k3                          |   ¦
% ---------------------------------------  x?--> z
%|     n1, k1 , thickness d             |
%|                                      |
% ---------------------------------------
%|     n2 , k2                          |
%|                                      |
%===================================================

b = linspace(k1,max(k2,k3),2000); % creates the \beta vector
% this loop calculates only the cases where the R.H.S. is real

kapa  = sqrt(k1^2 - b.^2); % variable of the Charecteristic Equation
gama  = sqrt(b.^2 - k2^2); % variabla of the Ch.Eq.
delta = sqrt(b.^2 - k3^2); % variabla of the Ch.Eq.

% L.H.S of equation (15)
f1 = tan(kapa*d);
% R.H.S of equation (15)
f2 = (kapa.*(gama + delta))./(kapa.^2 - delta.*gama);
%% Finding the TE modes
% Plots the two parts of the Ch.Eq. in order to find the specific points of
% intersection
f_1 = figure;
hold on
plot(b/k1,f1/max(f2),'Color','blue','Linewidth',1.5) % L.H.S.
plot(b/k1,f2/max(f2),'Color','red','Linewidth',1.5)  % R.H.S.
xlabel('\beta /\kappa_1')
ylabel('A.U.')
legend('L.H.S.','R.H.S.')
title('TE modes')
set(gca, 'XDir','reverse')
box on
ylim([-1 1]*1e-3)
f_1.Units = 'normalized';
f_1.OuterPosition = [0.1145    0.3081    0.3367    0.4756];

%%
% solver's options
options = optimset('Display','off');
% Characteristic equation
fun1 = @(t)tan(sqrt(k1^2 - t^2)*d) - (sqrt(k1^2 - t^2)*...
    (sqrt(t^2 - k2^2) + sqrt(t^2 - k3^2)))/(k1^2 - t^2 - sqrt(t^2 - k3^2)*sqrt(t^2 - k2^2));
neff = NaN*ones(1,TE_modes);
% Calculates the infinities of Tan(kappa*d)
for  n = 0:TE_modes-1
    k(n+1) = sqrt(k1^2 - ((pi/2/d)^2)*(2*n + 1)^2);
    neff(n+1) = fsolve(fun1,k(n+1)*(0.99999),options)/k0;
end
% Prints the results. Effective Refr. Index 
% and transmission constant
fprintf('\nThe results from the above figure are \n \n')
for i = 1:TE_modes
    fprintf('----> TE%d n_eff = %6.4f   beta = %6.2f (rad/um) <---- \n',...
        i-1,neff(i),neff(i)*k0*10^-6);
end
fprintf('\n')

%% Calculates N_eff of every guided mode
di  = (550:-5:1)*1e-9;         % span of thickness to be tested
xi = (n2^2-n3^2)/(n1^2-n2^2);  % APDWG "asymmetry" factor (0==SPDWG)
V  = k0*max(di)/2*sqrt( n1^2-n2^2 ); % V-number, same definition as for SPDWG
VasosTE = atan( sqrt(xi) )/2;  % V-number "asymmetry-offset" for TE-modes
TE_modes = ceil( (V-VasosTE)/(pi/2) ); % Number of modes expected
out = NaN*ones(TE_modes,length(di)); % initialization

k1_ = k1;
k1  = n1 *k0/k1_;       % wave-number in core
k2  = n2 *k0/k1_;       % wave-number in substrate
k3  = n3 *k0/k1_;       % wave-number in cladding layer

for u = 1:length(di)
    d = di(u);        % thickness tested
    V  = k0*d/2*sqrt( n1^2-n2^2 ); % V-number, same definition as for SPDWG
    TEmodes = ceil( (V-VasosTE)/(pi/2) ); % Number of modes expected
    % solver's options
    options = optimset('Display','off');
    % Characteristic equation to be solved for any d
    fun1 = @(t)tan(sqrt(k1^2 - t^2)*d*k1_) - (sqrt(k1^2 - t^2)*...
        (sqrt(t^2 - k2^2) + sqrt(t^2 - k3^2)))/(k1^2 - t^2 - sqrt(t^2 - k3^2)*sqrt(t^2 - k2^2));
    % Calculate the inifities of Tan(kappa * d)
    for  n = 0:TEmodes-1
        k(n+1) = sqrt(k1_^2 - ((pi/2/d)^2)*(2*n + 1)^2)/k1_;
    end
    % Solves the Ch.Eq. and calculates transmission constant \beta for
    % every guided TE mode. Also, a recursive method is used when is needed
    for n = 1:TEmodes
        ii = -2; % computational parameter
        out(n,u) = fsolve(fun1,k(n)*(0.99999),options);
        if (~isreal(out(n,u)) || out(n,u)<=k2 || out(n,u)>=k1)
            zz  = 1;% computational parameter
            flag = 0;% computational parameter
            while ((out(n,u)<=k2 || out(n,u)>=k1 || flag<1e-3) && ii<=-1)
                switch zz
                    case 1,
                        if (u==1)
                            out(n,u) = real(out(n,u));
                        elseif(isreal(out(n,u-1)))
                            out(n,u) = fsolve(fun1,out(n,u-1)*0.999,options);
                        else
                            out(n,u) = fsolve(fun1,k(n)*(1+0.1*10^ii),options);
                        end
                    case 2,
                        flag = out(n,u);
                        if(isreal(out(n,u)))
                            out(n,u)= (out(n,u) + fsolve(fun1,out(n,u)*(1+10*eps),options))/2;
                        else
                            out(n,u)= (out(n,u) + fsolve(fun1,k(n)*(1+0.1*10^ii),options))/2;
                        end
                        flag = abs(out(n,u) - flag);
                end
                zz = 2;% computational parameter
                ii= ii + 0.1;% computational parameter
            end
        end
    end
end
% Corrects the computational faults, if something went wrong.
% Ussually, this happens when the thickness step is too small
for n = 1:TE_modes
    for u = 2:length(di)-1
        if( ~isnan(out(n,u-1)*out(n,u+1)))
            if out(n,u)<out(n,u+1)
                mid = real((out(n,u-1) + out(n,u+1)))/2;
                out(n,u) = mid;
            end
        end
    end
end
n_eff = out*k1_/k0; % effective refractive index for every guided mode

%% Plot the results    N_eff vs Thickness d
f_2 = figure;
hold on; axis on; box on;
% TE_0
plot(di*10^9,real(n_eff(1,:)),'Color','blue','LineWidth',2);
% TE_1
plot(di*10^9,real(n_eff(2,:)),'Color','blue','LineWidth',2,'LineStyle','--');
% TE_2
plot(di*10^9,real(n_eff(3,:)),'Color','blue','Linewidth',2,'LineStyle','-.');

xlabel('Thickess (nm)')
ylabel('N_{eff}')
legend('TE_0','TE_1','TE_2')
ylim([n2*0.999 n1])
xlim([min(di) max(di)]*10^9)
f_2.Units = 'normalized';
f_2.OuterPosition = [0.4918    0.3069    0.3367    0.4756];
