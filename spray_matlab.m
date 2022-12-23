%                         ________  ________  ________  ________      ___    ___ ________                              %
%                        |\   ____\|\   __  \|\   __  \|\   __  \    |\  \  /  /|\   ___ \                             %
%                        \ \  \___|\ \  \|\  \ \  \|\  \ \  \|\  \   \ \  \/  / | \  \_|\ \                            %
%                         \ \_____  \ \   ____\ \   _  _\ \   __  \   \ \    / / \ \  \ \\ \                           %
%                          \|____|\  \ \  \___|\ \  \\  \\ \  \ \  \   \/  /  /   \ \  \_\\ \                          %
%                            ____\_\  \ \__\    \ \__\\ _\\ \__\ \__\__/  / /      \ \_______\                         %
%                           |\_________\|__|     \|__|\|__|\|__|\|__|\___/ /        \|_______|                         %
%                           \|_________|                            \|___|/                                            %
%                                                                                                                      %
%                                                                                                                      %
%                                            Author: Andrea Somma;                                                     % 
%                                            Politecnico of Milan 2021-2022                                            % 
%                                                                                                                      %     

clc; clear; close all

global P Gdry A B C rhoL Kpw nAvp h MWwater MWair dHev CpL CpG muG dp mp0 vg rho_milk fat kG Diff R Tg0

%% ---------------------------------------------------------------------------------------------------------------------
% Physical Properties
%% ---------------------------------------------------------------------------------------------------------------------
Tp0 = 303;								        % K
Tg0 = 403;										% K
P = 1;											% atm
CpL = 1;									    % kcal/kg/K
CpG = 0.25;										% kcal/kg/K
MWwater = 18;                                   % kg/kmol
MWair = 29;                                     % kg/kmol
R = 0.0821;                                     % m^3*atm/K/kmol
dHev = 540;										% kcal/kg
rhoL = 1000;					     			% kg/m^3
rhoG = P/R/Tg0*MWair;                           % kg/m^3
muG = 2.3*1e-5;                                 % kg/m/s
kG = 8*1.0e-6;									% kcal/m/s/K
Diff = 1.8*1e-5;                               	% m^2/s
rho_milk  =  440;                               % Kg/m3        

%% ---------------------------------------------------------------------------------------------------------------------
% Antoine equation
%% ---------------------------------------------------------------------------------------------------------------------
A = 18.3036;
B = 3816.44;
C = -46.13;

%% ---------------------------------------------------------------------------------------------------------------------
% Process stream
%% ---------------------------------------------------------------------------------------------------------------------
Qmilk = 1750/3600;                            	% kg/h -> kg/s
fat = 4.76/100;                           		% kg/kg
Qfat = Qmilk*fat;                     			% kg/s
Ql = Qmilk - Qfat;                				% kg/s
Win = Ql/Qfat;									% kg/kg
Wout = 0.005; 									% kg/kg
dp = 2*1e-4;     								% m
vp0 = 0.3;										% m/s
mp0 = rhoL * 3.14 / 6 * (dp^3);         		% kg
mdry = mp0 / (1 + Win);                   		% kg
mout = mdry * (1 + Wout);                  		% kg
nAvp = Qmilk / mp0;                            	% drops/s

%% ---------------------------------------------------------------------------------------------------------------------
% Gas stream
%% ---------------------------------------------------------------------------------------------------------------------
Gdry = 20;										% kg/s
D = 5.5;     									% m

%% ---------------------------------------------------------------------------------------------------------------------
% Momentum Balance
%% ---------------------------------------------------------------------------------------------------------------------
vg = Gdry / rhoG / 3.14 / D^2*4;                % m/s
vs = (rhoL-rhoG) * dp^2 * 9.81 / 18 / muG;      % m/s
vp = vg + vs;                                   % m/s
vs0 = vp0 - vg;                                 % m/s

%% ---------------------------------------------------------------------------------------------------------------------
% Mass and heat transfer
%% ---------------------------------------------------------------------------------------------------------------------
Re = rhoG * vs * dp / muG;                      % Reynolds
Pr = muG * CpG /kG; 						    % Prandtl
Sc =  muG / rhoG / Diff;      			    	% Schmidt
Nu = 2 + 0.4 * Re^0.5 * Pr^(1/3);          		% Nusselt
Sh = 2 + 0.4 * Re^0.5 * Sc^(1/3);              	% Sherwood
h = Nu * kG / dp;                          		% kcal/m^2/s/K
Kc = Sh * Diff / dp;                   			% m/s
Kpw = Kc / R / Tg0 * MWwater;      				% kg/m^2/s/atm

%% ---------------------------------------------------------------------------------------------------------------------
% Integration
%% ---------------------------------------------------------------------------------------------------------------------
tspan = linspace(0,5.5,1001);
y0 = [mp0 0 Tp0 Tg0 vs0 0];

[t,y] = ode15s(@solver,tspan,y0,odeset("MaxStep",1e-4,"Refine",5));

%% ---------------------------------------------------------------------------------------------------------------------
% Results
%% ---------------------------------------------------------------------------------------------------------------------
m = y(:,1);
time = t(find(m<2e-10,1));
lenght = y(find(m<2e-10,1),6);

%% ---------------------------------------------------------------------------------------------------------------------
% Plots
%% ---------------------------------------------------------------------------------------------------------------------
cc = winter(6);
close all

subplot(2,2,1)
plot(y(:,6), y(:,1),"LineWidth",1.4,"Color",cc(1,:))
hold on
plot([lenght lenght],[-1000 1000],"LineStyle","-.")
xlabel('Axial coordinate [m]','FontWeight','bold','FontSmoothing','on','FontSize',14)
ylabel('Particle mass [kg]','FontWeight','bold','FontSmoothing','on','FontSize',14)
ylim([min(y(:,1))/5,max(y(:,1))*1.1])
legend('Particle Mass','Usefull Lenght')
hold off

subplot(2,2,2)
plot(y(:,6),y(:,3),"LineWidth",1.4,"Color",cc(1,:))
hold on
plot(y(:,6),y(:,4),"LineWidth",1.4,"Color",cc(4,:))
plot([lenght lenght],[-1000 1000],"LineStyle","-.")
xlabel('Axial coordinate [m]','FontWeight','bold','FontSmoothing','on','FontSize',14)
ylabel('Temperature [K]','FontWeight','bold','FontSmoothing','on','FontSize',14)
legend('Particle','Air','Usefull Lenght')
ylim([min(y(:,3))/1.05,max(y(:,4))*1.05])
hold off

subplot(2,2,3)
plot(y(:,6),y(:,5),"LineWidth",1.4,"Color",cc(1,:))
hold on
plot(y(:,6), vg + y(:,5),"LineWidth",1.4,"Color",cc(4,:))
plot(y(:,6), vg * ones(1 * length(tspan)),"LineWidth",1.4,"Color",cc(6,:))
plot([lenght lenght],[-1000 1000],"LineStyle","-.","Color","r")
xlabel('Axial coordinate [m]','FontWeight','bold','FontSmoothing','on','FontSize',14)
ylabel('Velocity [m/s]','FontWeight','bold','FontSmoothing','on','FontSize',14)
legend('vs','vp','vg')
ylim([min(y(:,5))/5,max(vg + y(:,5))*1.1])
hold off

subplot(2,2,4)
plot(tspan,y(:,6),"LineWidth",1.4,"Color",cc(1,:))
hold on
plot([time time],[-1000 1000],"LineStyle","-.")
xlabel('Time [s]','FontWeight','bold','FontSmoothing','on','FontSize',14)
ylabel('Axial coordinate [m]','FontWeight','bold','FontSmoothing','on','FontSize',14)
legend('Axial Position','Usefull Lenght')
ylim([min(y(:,6)),max(y(:,6))*1.1])
hold off

%% ---------------------------------------------------------------------------------------------------------------------
% Function
%% ---------------------------------------------------------------------------------------------------------------------
function spraydryer = solver(t,x)

global P Gdry A B C rhoL Kpw nAvp h MWwater MWair dHev CpL CpG muG dp mp0 vg rho_milk fat kG Diff R
    
    mp = x(1);
    Gvap = x(2);
    Tp = x(3);
    Tg = x(4);
    vs = x(5);
    z = x(6);
    
    mmin = mp0 * fat;
    mp = max(mmin, mp);
    
    x_H2O = (mp - mp0*fat) / mp;
    
    rho_avg = x_H2O * rhoL + (1 - x_H2O) * rho_milk;
    V = mp / rho_avg;
    dp = (6 * V / pi)^(1/3);
    Dmin = 83.87e-6;
    dp = max(Dmin, dp);
    
    rhoG = P / R / Tg * MWair; 
    Re = abs(rhoG * vs * dp / muG);                     % Reynolds
    Pr = muG * CpG / kG; 								% Prandtl
    Sc = muG / rhoG / Diff;      						% Schmidt
    Nu = 2 + 0.4 * Re^0.5 * Pr^(1/3);                  	% Nusselt
    Sh = 2 + 0.4 * Re^0.5 * Sc^(1/3);                  	% Sherwood
    f = (1 + 0.14*Re^0.7);                              % Drag factor
    h = Nu * kG / dp;                             		% kcal/m^2/s/K
    Kc = Sh * Diff/ dp;                      			% m/s
    Kpw = Kc / R / Tg * MWwater;         			    % kg/m^2/s/atm
    
    Pw = MWair / MWwater * Gvap * P / Gdry;
    P0 = exp(A - B / (Tp + C)) / 760;                   % mmHg--->atm
    Sp = pi * dp^2;
    
    if mp > mmin + mmin/10000
        spraydryer(1) = Kpw * Sp * (Pw - P0);
    else
        spraydryer(1) = 0;
    end
    
    spraydryer(2) = -Kpw * Sp * (Pw-P0) * nAvp;
    spraydryer(3) = (h * Sp * (Tg - Tp) + spraydryer(1) * dHev) / mp / CpL;
    spraydryer(4) = -h * Sp * (Tg - Tp) * nAvp / (Gvap + Gdry) / CpG;
    spraydryer(5) = ((1 - rhoG / rho_avg) * 9.81 - 3 * f * vs * muG * 3.14 * dp /mp - vs * spraydryer(1) / mp);
    spraydryer(6) = (vs+vg);
    
    spraydryer = spraydryer';
end
