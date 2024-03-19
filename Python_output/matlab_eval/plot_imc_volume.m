% -----------------------------------
% Experimental velocity vs thickness curve
% Written by Erik Jacobsson 2023-02-07
% -----------------------------------

close all; clear; clc;

% Unit scale
mm_to_nm = 1e6;

% % Plot thickness vs time
% figure()
% plot(tvec_exp,dvec_exp*mm_to_micron,'LineWidth',2.0)
% xlabel('Time [hr]')
% ylabel('Thickness [\mu m]')
% title('Thickness vs. Time')
% 
% % Plot thickness velocity vs time
% figure()
% plot(tvec_exp,vvec_exp*mm_to_micron,'LineWidth',2.0)
% xlabel('Time [hr]')
% ylabel('Thickness velocity [\mu m/s]')
% title('Velocity vs. Time')

% % Plot thickness velocity vs thickness
% figure()
% plot(dvec_exp*mm_to_micron,vvec_exp*mm_to_micron,'LineWidth',2.0)
% xlabel('Thickness [\mu m]')
% ylabel('Thickness velocity [\mu m/hr]')
% title('Velocity vs. Thickness')
% axis([0 0.25 0 9e-3])


% % Plot thickness velocity vs thickness with rounded coefficient
% figure()
% dvec2 = linspace(0.001/mm_to_micron,0.22/mm_to_micron)';
% % vvec2 = 0.6*0.014^(1/0.6)*dvec2.^(-2/3)*1e-5;  % Observe!
% vvec2 = (4.88*1e-4)*dvec2.^(-2/3)*1e-5;  % Observe!
% vvec23 = -dvec2.^(-1)*(6.575e-19)*(1000)^2*3600*(42.3577-37.3853)/(42.3577-61.4786); % mol/mm3  ; % mol/mm3*dvec2.^(-1)*1e-5;  % Observe!
% treshhold = 5e-6;
% vvec2(vvec2>treshhold) = treshhold;
% vvec23(vvec23>10e-6) = 10e-6;
% plot(dvec2*mm_to_micron,vvec2*mm_to_micron,'LineWidth',2.0)
% hold on
% plot(dvec2*mm_to_micron,vvec23*mm_to_micron,'LineWidth',2.0)
% xlabel('Thickness [\mu m]')
% ylabel('Thickness velocity [\mu m/hr]')
% title('Velocity vs. Thickness')
% legend('d-2/3','d-1')
% % axis([0 0.25 0 9e-3])


%% Load simulation results

% % Matlab
% load(['/home/er7128ja/Nextcloud/Projekt/Level_set/Level_set_Matlab', ...
%      '/mat_files/level_set_IMCvol_2D4.mat'])
% tvec_4M = tvec; IMC_vol_per_area_4M = IMC_vol_per_area;

% MC1
load(['/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output/', ...
      'matlab_eval/mat_eval_files/IMC_volume_MC1.mat'])
tvec_MC1 = tvec; IMC_vol_per_area_MC1 = IMC_vol_per_area;


% MC2
load(['/home/er7128ja/Nextcloud/Projekt/Project_Code_2D/Python_output/', ...
      'matlab_eval/mat_eval_files/IMC_volume_MC2.mat'])
tvec_MC2 = tvec; IMC_vol_per_area_MC2 = IMC_vol_per_area;

%% Experimental (Brown)
tvec_exp_end = tvec(end);
tvec_exp_p = linspace(3,tvec_exp_end-40,6)';                % time [hr]
dvec_exp_p = 0.014*tvec_exp_p.^0.6*1e-3;        % thickness [mm]

tvec_exp = [3 8 21 38 67 95 120];
dvec_exp = [0.025 0.04 0.08 0.11 0.16 0.24 0.3]*1e-3;

% vvec_exp = 0.6*0.014*tvec_exp.^(-0.4)*1e-3; % velocity [mm/hr]

%% Plot IMC vol per area vs time for experimental and simulation results
figure()

% Simulation MC1
plot(tvec_MC1,IMC_vol_per_area_MC1*mm_to_nm,'k-','LineWidth',3.0,...
    'DisplayName','MC1')
hold on


% % Simulation MC2
% plot(tvec_MC2,IMC_vol_per_area_MC2*mm_to_nm,'-s','LineWidth',3.0,...
%     'Color',[0.4940 0.1840 0.5560], 'DisplayName','MC2',...
%     'MarkerSize',10,...
%     'MarkerFaceColor',[0.4940 0.1840 0.5560],...
%     'MarkerIndices',1:5:length(tvec_MC2))

% Experiment
plot(tvec_exp,dvec_exp*mm_to_nm,'^','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[1 .6 .6], 'DisplayName','Chason 2008')


plot(tvec_exp_p,dvec_exp_p*mm_to_nm,'>','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[1 .2 .2], 'DisplayName','Chason 2008 fitted')



%% Legend and axis

% Legend and axis
legend('Location','southeast')
legend boxoff 
axis([0,600,0,900])
% ylabel('$\mathbf{Cu_6Sn_5}$ \textbf{vol./area (nm)}','Interpreter','latex')
% xlabel('\textbf{t (h)}','Interpreter','latex')
xlabel('t (h)')
ylabel('Cu6Sn5volumeperarea')
ax = gca;
ax.FontSize = 14;
legend('Location','southeast')
legend boxoff 
 
% filename = '/home/er7128ja/Nextcloud/Projekt/Project_Code/articel_final/results/result_figs/imc_volume.eps';
% exportgraphics(gca,filename,'Resolution',300)

%% Save

% Create plot dir in single_study if not exist
plotfolder= '../single_study/plots';
if ~exist(plotfolder,"dir")
    mkdir(plotfolder)
    disp(['Created new directory: ', plotfolder])
end

filename = 'imc_volume.png';
saveas(gcf, fullfile(plotfolder,filename));
