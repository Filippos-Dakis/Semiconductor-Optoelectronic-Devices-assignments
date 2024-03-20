% Filippos Tzimkas-Dakis  UoC October 2020
% Homework 1 ------- Optoelectronics and Lasers
% This is a script to calculate the energy eigenvalues of a FINITE square
% well potential. 
% Also this script calculates the appropriate width of a FINITE square well
% potential for a fixed value of energy and mass.

% Solving the square welll potential 
%------ Even solutions -------------
eV2J  = 1.60218e-19;            % converts eVolts to Joules
J2eV  = 6.242e+18;              % converts Joules to eVolts
c     = 3*1e8;                  % speed of light (m/s)
Vo    = 0.3;                        % well's potential in (eV)
Vo    = Vo*eV2J;
me    = 9.1093837015*1e-31;         % electron mass in (Kg)
m     = 0.067*me;                       % effective mass of the electron 
                                    %  inside the potential in (Kg)
h_eV = 6.582119569*1e-16;           % h-bar in (eV*s)
h_J  = 1.054571817*1e-34;           % h-bar in (J*s)
% h_ev^2 /me 
L     = 10*1e-9;                    % width (m)

%% Firstly, we will find the solution graphically.
Uo = 2*m*Vo/h_J^2;             % defines a new normalized potential value
y  = sqrt(Uo)*L/2;             
z  = linspace(0.1,y,1e4);   % z = sqrt(2*me*E/h_J^2)*L/2
f1 = tan(z);                   % Function 1 of the trascedental equation
f2 = sqrt(y^2 - z.^2)./z;      % Function 2 of the trascedental equation

close all

figure(1)
plot(z,f1,'LineWidth',1.5,'Color','blue')
hold on
plot(z,f2,'LineWidth',1.5,'Color','red')
hold on
plot(z,z*0,'Color','black')
axis on
ylim([-0.1 20])
xlabel('Well`s width (nm)')
ylabel('Energy (eV)')

%% Pure Mathematical way for Even solutions
options = optimset('Display','off');
fun1 = @(t)tan(t) - sqrt(y^2 - t^2)./t;
out = fsolve(fun1,0.1,options);

E_e = (2*out/L)^2*(h_J^2 /2/ m)*J2eV  % even eneergy (eV)

%% Firstly, we will find the solution graphically.

% close all

figure(2)
plot(z,f1,'LineWidth',1.5,'Color','blue')
hold on
plot(z,- 1./f2,'LineWidth',1.5,'Color','red')
hold on
plot(z,z*0,'Color','black')
axis on
ylim([-20 0.1])
xlabel('Well`s width (nm)')
ylabel('Energy (eV)')

%% Pure Mathematical way for Odd solutions
options = optimset('Display','off');
fun1 = @(t)tan(t) + t/sqrt(y^2 - t^2);
out = fsolve(fun1,pi/2 + 0.1,options)

E_o = (2*out/L)^2*(h_J^2 /2/ m)*J2eV  % odd eneergy (eV)


