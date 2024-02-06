% --- Plot mean biaxial stress ---

close all; clc; clear; 


%% Load

% MC1
load(['/home/er7128ja/Nextcloud/Projekt/Project_Code/Python_output/', ...
      'matlab_eval/mat_eval_files/biax_stressMC1.mat'])
tvec_MC1 = tvec; biax_Sn_avg_MC1 = biax_Sn_avg;


% MC2
load(['/home/er7128ja/Nextcloud/Projekt/Project_Code/Python_output/', ...
      'matlab_eval/mat_eval_files/biax_stressMC2.mat'])
tvec_MC2 = tvec; biax_Sn_avg_MC2 = biax_Sn_avg;

%% Plot


% Plot mean biaxial stress vs time
f1 = figure();

% MC1
plot(tvec_MC1, biax_Sn_avg_MC1,'k-','LineWidth',3.0, 'DisplayName','MC1')
hold on

% MC2
plot(tvec_MC2,biax_Sn_avg_MC2,'-s','LineWidth',3.0,...
    'Color',[0.4940 0.1840 0.5560], 'DisplayName','MC2',...
    'MarkerSize',10,...
    'MarkerFaceColor',[0.4940 0.1840 0.5560],...
    'MarkerIndices',1:5:length(tvec_MC2))


% Exp Lee
tlee = [ 2  40   50  60  150  170 220 280 320 450 520];
slee = [-1 -3.5 -4  -4  -7    -7  -7  -7  -8  -7 -7.5];
plot(tlee,slee,'>','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[1 .2 .2], 'DisplayName','Lee 2002')

% Exp brown
tbrown = [8 20 38 65 90 120];
sbrown = [-12 -8 -11 -14 -13 -14];
plot(tbrown,sbrown,'^','MarkerSize',10,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor',[1 .6 .6], 'DisplayName','Chason 2008')

% Labels and axis
% xlabel('\textbf{t (h)}','Interpreter','latex')
% ylabel('\textbf{Mean biaxial stress (MPa)}','Interpreter','latex');
xlabel('t (h)')
ylabel('meanbiaxialstressMPaAt')
% axis([0,tvec(k-1),-15,7])
axis([0,600,-15,1])
ax = gca;
ax.FontSize = 14;
legend('Location','southeast')
legend boxoff 


% Save results in plot dir in single_study if not exist
filename = '/home/er7128ja/Nextcloud/Projekt/Project_Code/articel_final/results/result_figs/biax_avg.eps';
exportgraphics(gca,filename,'Resolution',300)



%% Functions


function save_plot(fh,filename)
plotfolder= '../single_study/plots';
if ~exist(plotfolder,"dir")
    mkdir(plotfolder)
    disp(['Created new directory: ', plotfolder])
end
exportgraphics(fh,fullfile(plotfolder,filename),'Resolution',150)
end