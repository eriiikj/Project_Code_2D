% -----------------------------------
% Fitted Gibbs free energy curves
% Written by Erik Jacobsson 2023-03-29
% -----------------------------------

close all; clear; clc;

% Article colors, red, orange, blue, green
RGB = [[255,58,58]; [255,165,41]; [8,81,156]; [255,58,58]]/255;

% Molar mass [g/mol]
M_Cu  = 63.55;
M_Sn  = 118.71;
x_Sn  = 6/11;
M_IMC = x_Sn*M_Cu + (1-x_Sn)*M_Sn;

% Molar volumes [mm3/mol]
V_Cu  = 7.09*1e3;
V_IMC = 10.7*1e3;
V_Sn  = 16.3*1e3;


% Mean compositions (xbar)
xh_cu  = 0.10569;
xh_imc = 0.41753;
xh_sn  = 0.99941;

% Free energy parameters (A, B, C)
cu_params  = [1.0133e5 -2.1146e4 -1.2842e4 xh_cu]';
imc_params = [4e5      -6.9892e3 -1.9185e4 xh_imc]';
sn_params  = [4.2059e6  7.1680e3 -1.5265e4 xh_sn]';


% Composition array
xvec = linspace(0,1,5000)';


% Fitted parabola of form Gi = Ai/2(xi-xih)^2 + Bi(xi-xih) + Ci
G_cu  = get_G(xvec, cu_params);
G_imc = get_G(xvec, imc_params);
G_sn  = get_G(xvec, sn_params);

% Fitted form of dGi/dx
dGdx_cu  = get_dGdx(xvec, cu_params);
dGdx_imc = get_dGdx(xvec, imc_params);
dGdx_sn  = get_dGdx(xvec, sn_params);

tol     = 1e-7;
% --- Cu/IMC ---
x_imccu = [0.1;0.4];
x_imccu = rhaps(x_imccu,1e-4,imc_params,cu_params)

% Tangent array for cu/imc at x0(1)
xvec1 = linspace(0.03,0.44,5000)';
dGvec_cuimc1 = get_dGvec(x_imccu(1),xvec1,imc_params);
P_cuimc1     = [x_imccu(1) get_G(x_imccu(1), imc_params)];

% Tangent array for cu/imc at x0(2)
xvec2 = linspace(0.03,0.44,5000)';
dGvec_cuimc2 = get_dGvec(x_imccu(2),xvec2,cu_params);
P_cuimc2     = [x_imccu(2) get_G(x_imccu(2), cu_params)];

% --- Sn/IMC ---
x_snimc = [0.8; 0.6];
x_snimc = rhaps(x_snimc,1e-4,sn_params,imc_params)

% Tangent array for imc/sn at x0(1)
xvec3 = linspace(0.38,1,5000)';
dGvec_snimc1 = get_dGvec(x_snimc(1),xvec3,sn_params);
P_imcsn1     = [x_snimc(1) get_G(x_snimc(1), sn_params)];

% Tangent array for imc/sn at x0(2)
xvec4 = linspace(0.38,1,5000)';
dGvec_snimc2 = get_dGvec(x_snimc(2),xvec4,imc_params);
P_imcsn2     = [x_snimc(2) get_G(x_snimc(2), imc_params)];


% --- Cu/Sn ---
x_sncu = [0.8; 0.2];
x_sncu = rhaps(x_sncu,1e-4,sn_params,cu_params)

% Tangent array for imc/sn at x0(1)
xvec5 = linspace(0.16,1,5000)';
dGvec_sncu1 = get_dGvec(x_sncu(1),xvec5,sn_params);
P_sncu1     = [x_sncu(1) get_G(x_sncu(1), sn_params)];

% Tangent array for imc/sn at x0(2)
xvec5 = linspace(0.16,1,5000)';
dGvec_sncu2 = get_dGvec(x_sncu(2),xvec5,cu_params);
P_sncu2     = [x_sncu(2) get_G(x_sncu(2), cu_params)];

% --- Plot ---
% Plot free energies
figure()
plot(xvec,G_cu,'color',RGB(1,:),'LineStyle','-','LineWidth',3,'DisplayName','G1')
hold on
plot(xvec,G_imc,'color',RGB(2,:),'LineStyle','--','LineWidth',3,'DisplayName','G2')
plot(xvec,G_sn,'color',RGB(3,:),'LineStyle',':','LineWidth',3,'DisplayName','G3')


% Plot tangent line for cu/imc
plot(xvec1,dGvec_cuimc1,'k-','LineWidth',1,'DisplayName','$dG_{Cu}/dx$', ...
                     'HandleVisibility','off')
plot(xvec2,dGvec_cuimc2,'k-','LineWidth',1,'DisplayName','$dG_{Cu}/dx$', ...
                     'HandleVisibility','off')
plot(P_cuimc1(1),P_cuimc1(2),'ko','HandleVisibility','off', 'MarkerSize',8)
plot(P_cuimc2(1),P_cuimc2(2),'ko','HandleVisibility','off', 'MarkerSize',8)

% Plot tangent line for imc/sn
plot(xvec3,dGvec_snimc1,'k-','LineWidth',1,'DisplayName','$dG_{Cu}/dx$', ...
                     'HandleVisibility','off')
plot(xvec4,dGvec_snimc2,'k-','LineWidth',1,'DisplayName','$dG_{Cu}/dx$', ...
                     'HandleVisibility','off')
plot(P_imcsn1(1),P_imcsn1(2),'ko','HandleVisibility','off', 'MarkerSize',8)
plot(P_imcsn2(1),P_imcsn2(2),'ko','HandleVisibility','off', 'MarkerSize',8)

% Plot tangent line for cu/imc
plot(xvec5,dGvec_sncu1,'k-','LineWidth',1,'DisplayName','$dG_{Cu}/dx$', ...
                     'HandleVisibility','off')
plot(xvec5,dGvec_sncu2,'k-','LineWidth',1,'DisplayName','$dG_{Cu}/dx$', ...
                     'HandleVisibility','off')
plot(P_sncu1(1),P_sncu1(2),'ko','HandleVisibility','off', 'MarkerSize',8)
plot(P_sncu2(1),P_sncu2(2),'ko','HandleVisibility','off', 'MarkerSize',8)


% Legend and axis
legend('Location','best')
legend boxoff 
axis([0 1 -20000 -5000])
% xlabel('$\textbf{Molar fraction of Sn (-)}$','Interpreter','latex')
xlabel('Molar fraction of Sn (-)')
% ylabel('\textbf{Gibbs energy (Jmol}$\boldmath{^{-1}}\textbf{)}$','Interpreter','latex')
ylabel('Gibbs energy (Jmol1)as')
ax = gca;
ax.FontSize = 14;
ax.YAxis.Exponent = 0;
ytickformat('%.0f')
set(gcf,'renderer','Painters')
saveas(gca, '/home/er7128ja/Nextcloud/Projekt/Project_Code/articel_final/figs/gibbs_energies','epsc')


% % Determine diffusion potentials (mu = A(x-xh) + B)
% mu_cuimc = cu_params(1)*(x_imccu(2) - cu_params(4)) + cu_params(2)
% mu_imccu = imc_params(1)*(x_imccu(1) - imc_params(4)) + imc_params(2)
% 
% mu_snimc = sn_params(1) *(x_snimc(1) - sn_params(4)) + sn_params(2)
% mu_imcsn = imc_params(1)*(x_snimc(2) - imc_params(4)) + imc_params(2)
% 
% 
% % Corresponding concentrations
% xcu_imc = (mu_cuimc - cu_params(2))/cu_params(1) + cu_params(4);
% ximc_cu = (mu_imccu - imc_params(2))/imc_params(1) + imc_params(4);
% ccu_imc = xcu_imc/V_Cu % [mol/mm3]
% cimc_cu = ximc_cu/V_IMC % [mol/mm3]
% 
% 
% % Corresponding concentrations
% xsn_imc = (mu_snimc - sn_params(2))/sn_params(1) + sn_params(4);
% ximc_sn = (mu_imcsn - imc_params(2))/imc_params(1) + imc_params(4);
% csn_imc = xsn_imc/V_Sn  % [mol/mm3]
% cimc_sn = ximc_sn/V_IMC % [mol/mm3]





%% Functions

% --- Gibbs free energy---
function G = get_G(x, fparams)
xh = fparams(4);
G  = fparams(1)/2*(x-xh).^2 + fparams(2)*(x-xh) + fparams(3);
end

% --- Derivative of Gibbs free energy ---
function dGdx = get_dGdx(x, fparams)
xh   = fparams(4);
dGdx = fparams(1)*(x-xh) + fparams(2);
end


function Kt = get_K(x,fparams1,fparams2)
x1  = x(1); x2 = x(2);
xh1 = fparams1(4); xh2 = fparams2(4);
A1  = fparams1(1); B1 = fparams1(2); 
A2  = fparams2(1); B2 = fparams2(2);
k1  = get_dGdx(x1, fparams1);
k2  = get_dGdx(x2, fparams2);

% Tangent
Kt  = [A1, -A2; (k1-(A1*(2*x1-xh1)+B1)) (k2-(A2*(2*x2-xh2)+B2))];
end


function r = get_r(x,fparams1,fparams2)
x1 = x(1); x2 = x(2);
k1 = get_dGdx(x1, fparams1);
k2 = get_dGdx(x2, fparams2);
m1 = get_G(x1, fparams1) - k1*x1;
m2 = get_G(x2, fparams2) - k2*x2;

% Residual
r = [k1 - k2; m1 - m2];
end

function x = rhaps(x0,tol,fparams1,fparams2)
res  = 1e6;
step = 1;
while (res>tol && step<=300)
r    = get_r(x0,fparams1,fparams2);    
Kt   = get_K(x0,fparams1,fparams2);
dx   = Kt\(-r);
x0   = x0 + dx;
res  = norm(r);
% disp(['res: ',num2str(res)])
step = step + 1;
end
x = x0;
end


function dGvec = get_dGvec(x,xvec,fparams)
% Tangent array
k     = get_dGdx(x,fparams);
m     = get_G(x,fparams) - k*x;
dGvec = k*xvec+m;
end