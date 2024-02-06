% -----------------------------------
% Fitted Gibbs free energy curves
% Written by Erik Jacobsson 2023-03-29
% -----------------------------------

close all; clear; clc;

% Sn composition
xsn       = linspace(0,1,100)';

% Fitted parabola of form Gi = Ai/2(xi-xih)^2 + Bi(xi-xih) + Ci
gibbs_cu  = get_gibbs_cu(xsn);
gibbs_imc = get_gibbs_imc(xsn);
gibbs_sn  = get_gibbs_sn(xsn);

% Tangents dGi/dx(x1) = Ai(xi-xih) + Bi
x = 0.10569;
k  = get_dgibbs_cu(x);
m  = get_gibbs_cu(x)-k*x;
dgibbs_cux1 = k*xsn+m;

% Tangents dGi/dx(x) = Ai(xi-xih) + Bi
x = 0.45290;
k  = get_dgibbs_imc(x);
m  = get_gibbs_imc(x)-k*x;
dgibbs_imcx2 = k*xsn+m;

% Tangents dGi/dx(x7) = Ai(xi-xih) + Bi
x5 = 0.312;
k  = get_dgibbs_cu(x5);
m  = get_gibbs_cu(x5)-k*x5;
dgibbs_cux7 = k*xsn+m;

% Points
x1 = 0.10569; x2 = 0.38210; x3=0.45290; x4 = 0.99941;
P1 = [x1 get_gibbs_cu(x1)];
P2 = [x2 get_gibbs_imc(x2)];
P3 = [x3 get_gibbs_imc(x3)];
P4 = [x4 get_gibbs_sn(x4)];

P5 = [x5 get_gibbs_cu(x5)];

% Plot free energies
plot(xsn,gibbs_cu,'-','LineWidth',2,'DisplayName','$$G_{Cu}$$')
hold on
plot(xsn,gibbs_imc,'-','LineWidth',2,'DisplayName','$$G_{IMC}$$')
plot(xsn,gibbs_sn,'-','LineWidth',2,'DisplayName','$$G_{Sn}$$')

% Plot tangent lines at equilibrium compositions
plot(xsn,dgibbs_cux1,'k-','LineWidth',1,'DisplayName','$$dG_{Cu}/dx$$', ...
                     'HandleVisibility','off')
plot(xsn,dgibbs_imcx2,'k-','LineWidth',1,'DisplayName','$$dG_{Cu}/dx$$',...
                      'HandleVisibility','off')

plot(xsn,dgibbs_cux7,'k-','LineWidth',1,'DisplayName','$$dG_{Cu}/dx$$',...
                      'HandleVisibility','off')

% Mark equilibrium compositions
plot(P1(1),P1(2),'ko','HandleVisibility','off')
plot(P2(1),P2(2),'ko','HandleVisibility','off')
plot(P3(1),P3(2),'ko','HandleVisibility','off')
plot(P4(1),P4(2),'ko','HandleVisibility','off')

plot(P5(1),P5(2),'ko','HandleVisibility','off')

% Legend and axis
legend('Location','southeast','Interpreter','latex','FontSize',12)
axis([0 1 -20000 -5000])
xlabel('Molar fraction of Sn')
saveas(gcf,'mat_eval_files/gibbs_energies.png')




%% Functions

% --- Cu Gibbs energy---
function gibbs_cu = get_gibbs_cu(x)
xh_cu     = 0.10569;
cu_params = [1.0133e5 -2.1146e4 -1.2842e4]';
gibbs_cu  = cu_params(1)/2*(x-xh_cu).^2 + cu_params(2)*(x-xh_cu) + ...
            cu_params(3);
end

function dgibbs_cu = get_dgibbs_cu(x)
xh_cu     = 0.10569;
cu_params = [1.0133e5 -2.1146e4 -1.2842e4]';
dgibbs_cu = cu_params(1)*(x-xh_cu) + cu_params(2);
end


% --- IMC Gibbs energy---
function gibbs_imc = get_gibbs_imc(x)
xh_imc     = 0.41753;
imc_params = [4e5 -6.9892e3 -1.9185e4]';
gibbs_imc  = imc_params(1)/2*(x-xh_imc).^2 + imc_params(2)*(x-xh_imc) + ...
             imc_params(3);
end

function dgibbs_imc = get_dgibbs_imc(x)
xh_imc     = 0.41753;
imc_params = [4e5 -6.9892e3 -1.9185e4]';
dgibbs_imc = imc_params(1)*(x-xh_imc) + imc_params(2);
end


% --- Sn Gibbs energy---
function gibbs_sn = get_gibbs_sn(x)
xh_sn     = 0.99941;
sn_params = [4.2059e6 7.1680e3 -1.5265e4]';
gibbs_sn  = sn_params(1)/2*(x-xh_sn).^2 + sn_params(2)*(x-xh_sn) + ...
            sn_params(3);
end

function dgibbs_sn = get_dgibbs_sn(x)
xh_sn     = 0.99941;
sn_params = [4.2059e6 7.1680e3 -1.5265e4]';
dgibbs_sn = sn_params(1)*(x-xh_sn) + sn_params(2);
end


