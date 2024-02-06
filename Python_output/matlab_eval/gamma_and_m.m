% Calculating values of the grain boundary energy gamma and mobility m

close all; clear; clc;
%% Previous values

% Mobility
mp_mm = 0.1;                   % mm^3/(Nh)
mp_m  = mp_mm*(1e-3)^3/(3600); % m^3/(Ns)

% Grain boundary energy
gammap_mm = 0.625*1e-6*(1e-4);  % J/mm^2
gammap_m = gammap_mm/(1e-3)^2;   % J/m^2


% Factor gamma*m
fac_c_prev   = gammap_mm*mp_mm
% fac_rho_prev = mp_mm


%% New model

disp('-------------')

% Mobility
beta     = 1;
deltaDgb = 7*1e-11;
b        = 0.256*1e-9;
Vm       = 10.4*(1e-2)^3;
R        = 8.31446261815324;
T        = 25 + 273.15;
Qm       = 104e3;
m0       = beta*deltaDgb*Vm/(b^2*R*T);
m        = m0*exp(-Qm/(R*T))
mgoal    = 1e-18

% Grain boundary energy
gamma = 0.625;

% Corresponding units in mm and hr
m_mm     = m*(1e3)^3/(1/3600);
gamma_mm = gamma/(1e3)^2;

fac_c   = gamma_mm*m_mm
% fac_rho = m_mm


