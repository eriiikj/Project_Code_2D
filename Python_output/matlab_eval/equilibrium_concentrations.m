% --- Computing equilibrium concentrations ---

clear; close all; clc;


% Mean compositions (xbar)
xh_cu  = 0.10569;
xh_imc = 0.41753;
xh_sn  = 0.99941;

% Free energy parameters (A, B, C)
cu_params  = [1.0133e5 -2.1146e4 -1.2842e4 xh_cu]';
imc_params = [4e5      -6.9892e3 -1.9185e4 xh_imc]';
sn_params  = [4.2059e6  7.1680e3 -1.5265e4 xh_sn]';

% Equilibrium compositions of Sn
x_cu_imc = 0.10571;
x_imc_cu = 0.38214;
x_imc_sn = 0.45292;
x_sn_imc = 0.99941;

% Molar mass [g/mol]
M_Cu  = 63.55;
M_Sn  = 118.71;
x_Sn  = 6/11;
M_IMC = x_Sn*M_Cu + (1-x_Sn)*M_Sn;

% Density [g/cm3]
rho_Cu  = 8.96;
rho_IMC = 8.28;
rho_Sn  = 7.30;

% Molar volumes [mm3/mol]
V_Cu  = 7.09*1e3;
V_IMC = 10.7*1e3;
V_Sn  = 16.3*1e3;

c_cu_imc = x_cu_imc/V_Cu;
c_imc_cu = x_imc_cu/V_IMC;
c_imc_sn = x_imc_sn/V_IMC;
c_sn_imc = x_sn_imc/V_Sn;


% Determine diffusion potentials (mu = A(x-xh) + B)
mu_cu_imc = cu_params(1)*(x_cu_imc  - cu_params(4))  + cu_params(2);
mu_imc_cu = imc_params(1)*(x_imc_cu - imc_params(4)) + imc_params(2);
mu_imc_sn = imc_params(1)*(x_imc_sn - imc_params(4)) + imc_params(2);
mu_sn_imc = sn_params(1) *(x_sn_imc - sn_params(4))  + sn_params(2);


% Diffusion coefficients (mm2/h)
DCu  = 2.877*10^-36*(1d6*3600);
DIMC = 6.575*10^-19*(1d6*3600);
DSn  = 2.452*10^-17*(1d6*3600);

% Mobilities
MCu  = DCu/(V_Cu*cu_params(1));
MIMC = DIMC/(V_IMC*imc_params(1));
MSn  = DSn/(V_Sn*sn_params(1));

% Observe equality
k1 = DIMC*(c_imc_sn - c_imc_cu);
k2 = MIMC*(mu_imc_sn-mu_imc_cu);