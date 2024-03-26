% -----------------------------------
% Extracting and plotting biax stress 
% from .mat files generated by whisker program.
% Written by Erik Jacobsson 2022-01-18
% -----------------------------------

close all; clc; clear

%% Import data

% Load level set init
[coord,enod,nodel,ex,ey,gpx,gpy,nelm,nnod,edof_1D,axisbc,IMC_steps,...
    ngrains] = load_level_set_init();

%% Biax stress in Sn
IMC_steps   = 180;
step_size   = 10;
niterations = floor((IMC_steps - 1) / step_size);
biax_Sn_avg = zeros(niterations,1);
tvec        = zeros(niterations,1);

k = 1;
for i_IMC=1:step_size:IMC_steps

    disp('--------')
    disp(['IMC step ', num2str(i_IMC)])

    % Load time
    tvec(k) = get_time(i_IMC);

    % Load stress state
    biax_Sn_avg(k) = load_biaxsnavg(i_IMC);  

    % Update counter
    k = k + 1;
end


figure(1)
plot(tvec, biax_Sn_avg,'k-','LineWidth',3.0, 'DisplayName','MC1')
hold on

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
xlabel('\textbf{t (h)}','Interpreter','latex')
ylabel('\textbf{Mean biaxial stress (MPa)}','Interpreter','latex');
% axis([0,tvec(k-1),-15,7])
axis([0,600,-15,1])
ax = gca;
ax.FontSize = 12;
legend('Location','southeast')
legend boxoff 


% Save mean biaxial stress
save('mat_eval_files/biax_stressMC1.mat','biax_Sn_avg','tvec')

%% Save plot

% Create plot dir in single_study if not exist
plotfolder= '../single_study/plots';
if ~exist(plotfolder,"dir")
    mkdir(plotfolder)
    disp(['Created new directory: ', plotfolder])
end

filename = 'biax_Sn_avg.png';
saveas(gcf, fullfile(plotfolder,filename));

%% Functions

function [coord,enod,nodel,ex,ey,gpx,gpy,nelm,nnod,edof_1D,axisbc,...
    IMC_steps,ngrains] = load_level_set_init()
% --- Function for loading initial level set data ---

% Filename
s             = what('../single_study/mat_files');
file_location = s.path;
filename      = 'level_setinit.mat';
filename      = fullfile(file_location,filename);

% Load arrays
load(filename,'coord')
load(filename,'enod')
load(filename,'nodel')
load(filename,'ex')
load(filename,'ey')
load(filename,'gpx')
load(filename,'gpy')
load(filename,'IMC_steps')
load(filename,'ngrains')

% Rearrange arrays
coord            = coord';
enod             = double(enod)';
nodel            = double(nodel);
ex               = ex';
ey               = ey'; 
gpx              = gpx';
gpy              = gpy';
gpx(:,[1 2 3 4]) = gpx(:,[1 2 4 3]);
gpy(:,[1 2 3 4]) = gpy(:,[1 2 4 3]);
IMC_steps        = double(IMC_steps);
ngrains          = double(ngrains);

% Extract other properties
nelm             = size(enod,1);
nnod             = max(enod(:));
edof_1D          = zeros(nelm,nodel+1);
edof_1D(:,1)     = 1:nelm;
edof_1D(:,2:end) = enod;
axisbc           = [min(ex(:)) max(ex(:)) min(ey(:)) max(ey(:))]';
end


function biax_Sn_avg = load_biaxsnavg(i_IMC)
s = what('../single_study/mat_files');
file_location = s.path;
filename = ['stress_',num2str(i_IMC),'.mat'];
filename = fullfile(file_location,filename);
load(filename,'biax_Sn_avg')
end

function time = get_time(itot)
s = what('../single_study/mat_files');
file_location = s.path;
filename = ['level_set_',num2str(itot),'.mat'];
filename = fullfile(file_location,filename);
load(filename,'time')
end