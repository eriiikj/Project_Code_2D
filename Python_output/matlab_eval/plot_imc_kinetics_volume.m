% --- Plotting IMC volume ---
clc; clear; close all;


%% Load
% Unit conversion
mm_to_um  = 1e3;
mm_to_nm  = 1e6;

% Load 

% Alpha=1d4
load(['/home/er7128ja/Nextcloud/Projekt/Project_Code/Python_output/', ...
      'matlab_eval/mat_eval_files/IMC_volume_alpha1d4.mat'])
tvec_1d4 = tvec; IMC_vol_per_area_1d4 = IMC_vol_per_area;

% Alpha=2d4
load(['/home/er7128ja/Nextcloud/Projekt/Project_Code/Python_output/', ...
      'matlab_eval/mat_eval_files/IMC_volume_alpha2d4.mat'])
tvec_2d4 = tvec; IMC_vol_per_area_2d4 = IMC_vol_per_area;

% Alpha=5d4
load(['/home/er7128ja/Nextcloud/Projekt/Project_Code/Python_output/', ...
      'matlab_eval/mat_eval_files/IMC_volume_alpha5d4.mat'])
tvec_5d4 = tvec; IMC_vol_per_area_5d4 = IMC_vol_per_area;

% Load
load('mat_eval_files/IMC_volume.mat')
h0 = 0;%IMC_vol_per_area(1);

%% --- Plot IMC_vol_per_area vs time ---

% Article colors, red, orange, blue, green
RGB = [[255,58,58]; [255,165,41]; [8,81,156]; [255,58,58]]/255;

% --- Plot IMC_vol_per_area vs time ---
f1 = figure();

% Plot analytical solution
tvec = linspace(0,tvec_1d4(end));
[t_al,h_al,t0] = get_analytical(tvec, IMC_vol_per_area);
plot(t_al,h_al*mm_to_nm,'k-','LineWidth',3.0,...
    'DisplayName','Analytical sol.')
hold on

% Plot simulation data
plot(tvec_1d4,IMC_vol_per_area_1d4*mm_to_nm,'color',RGB(1,:),'LineStyle',...
    '--','LineWidth',3.0,'DisplayName','$\alpha=1\times 10^4$')
plot(tvec_2d4,IMC_vol_per_area_2d4*mm_to_nm,'color',RGB(2,:),'LineStyle',...
    ':','LineWidth',3.0,'DisplayName','$\alpha=2\times 10^4$')
plot(tvec_5d4,IMC_vol_per_area_5d4*mm_to_nm,'color',RGB(3,:),'LineStyle',...
    '-.','LineWidth',3.0,'DisplayName','$\alpha=5\times 10^4$')

% Label
% xlabel('\textbf{t (h)}','Interpreter','latex')
xlabel('t(h)')
% ylabel('\textbf{Cu$$_6$$Sn$$_5$$ vol./area (nm)}','Interpreter','latex')
ylabel('Cu6Sn5volumearea')
legend('Location','southeast','Box','off','Interpreter','latex','FontSize',16)
axis([0, 600, 0, 300])
ax = gca;
ax.FontSize = 16;
yticks([0 100 200 300])

% Save results in plot dir in single_study if not exist
plotfolder= '../single_study/plots';
if ~exist(plotfolder,"dir")
    mkdir(plotfolder)
    disp(['Created new directory: ', plotfolder])
end
% set(gcf,'renderer','Painters')
% saveas(gca, '/home/er7128ja/Nextcloud/Projekt/Project_Code/articel_final/results/result_figs/interface_kin1','epsc')
filename = '/home/er7128ja/Nextcloud/Projekt/Project_Code/articel_final/results/result_figs/interface_kin1.eps';
exportgraphics(gca,filename,'Resolution',300)

% %% --- Plot IMC_vol_per_area vs sqrt(time) ---
f2 = figure();

% Plot analytical solution
t_al2 = t_al ;
plot(sqrt(t_al2),h_al*mm_to_nm,'k-','LineWidth',3.0,...
    'DisplayName','Analytical sol.')
hold on

% Plot simulation data
tvec2 = tvec_1d4;
plot(sqrt(tvec2),IMC_vol_per_area_1d4*mm_to_nm,'color',RGB(1,:),...
    'LineStyle','--','LineWidth',3.0,'DisplayName','$$\alpha=1\times 10^4$$')

% Label
% xlabel('\textbf{$$\mathbf{\sqrt{t}}$$ (h)}','Interpreter','latex')
xlabel('sqrtt(h)')
% ylabel('\textbf{Cu$$_6$$Sn$$_5$$ vol./area (nm)}','Interpreter','latex')
ylabel('Cu6Sn5volumearee')
legend('Location','southeast','Box','off','Interpreter','latex','FontSize',16)
% set(gca,'ytick',[])
set(gca,'yticklabel',[])
axis([0, 25, 0, 300])
ax = gca;
ax.FontSize = 16;
yticks([0 100 200 300])

% Save
% set(gcf,'renderer','Painters')
% saveas(gca, '/home/er7128ja/Nextcloud/Projekt/Project_Code/articel_final/results/result_figs/interface_kin2','epsc')
filename = '/home/er7128ja/Nextcloud/Projekt/Project_Code/articel_final/results/result_figs/interface_kin2.eps';
exportgraphics(gca,filename,'Resolution',300)


%% Functions
function [t_al,h_al,t0] = get_analytical(tvec, IMC_vol_per_area)

% Initial area
h0 = IMC_vol_per_area(1); % [mm]

% --- Analytical solution ---
t_al_end = tvec(end);
t_al     = linspace(0,t_al_end,2000)'; % [hr]

% Concentrations
cimc_cu = 3.5714e-05;
cimc_sn = 4.2329e-05;
csn_imc = 6.1313e-05;

% D
% DIMC       = 1d0*(6.575d-19)*(1d3)^2d0*3600d0; % [mm/h]
DIMC       = (3.575d-20)*(1d3)^2d0*3600d0;
k1         = 1/(csn_imc - cimc_sn);
k2         = DIMC*(cimc_sn - cimc_cu);
kd         = k1*k2;

% M
VmIMC      = 10.7*1e3;
AIMC       = 4e5;
MIMC       = DIMC/(AIMC*VmIMC);

muIMC_low  = -2.1144e+04;
muIMC_high = 7.1673e+03;
ka         = 1/(csn_imc - cimc_sn);
kb         = MIMC*(muIMC_high - muIMC_low);
km         = ka*kb;

% Height
k     = kd;
C     = 0;%h0^2;
h_al  = sqrt(2*k*t_al + C);


% Time axis
t0 = h0^2/(2*k);
end


function [tvec_exp,dvec_exp] =  plot_experimental(tvec)


% Experimental solution (Brown)
tvec_exp = linspace(0,tvec(end),10)';   % time [hr]
dvec_exp = 0.014*tvec_exp.^0.6*1e-3; % thickness [mm]

% % Plot experimental data
% plot(tvec_exp,(a1+dvec_exp)*mm_to_micron,'r-','LineWidth',2.0,...
%     'DisplayName','Experiment (Brown)')
% hold on

end

function save_plot(fh,filename)
plotfolder= '../single_study/plots';
if ~exist(plotfolder,"dir")
    mkdir(plotfolder)
    disp(['Created new directory: ', plotfolder])
end
exportgraphics(fh,fullfile(plotfolder,filename),'Resolution',150)
end