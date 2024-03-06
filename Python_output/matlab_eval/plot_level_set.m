% -----------------------------------
% Extracting and plotting level set properties  
% .mat files generated by Fortran program.
% Written by Erik Jacobsson 2023-09-05
% -----------------------------------

clear; close all; clc;

% Add path
addpath(genpath("Calfem/"))

% Load level set init
[coord,enod,nodel,ex,ey,gpx,gpy,nelm,nnod,edof_1D,axisbc,IMC_steps,...
    ngrains] = load_level_set_init();

% Length
L = axisbc(2)-axisbc(1);

% Get colors
my_colors = get_my_colors();

% Article colors, red, orange, blue, green
my_colors2 = [[255,58,58]; [255,165,41]; [8,81,156]; [136,86,167]
              [44,162,95]]/255;

% Unit conversion
mm_to_um  = 1e3;

% IMC steps
IMC_steps = 4;


%% Level set evolution
step_size   = 1;
niterations = floor((IMC_steps - 1) / step_size);

% Storing all times
tvec             = zeros(niterations,1);
IMC_vol_per_area = zeros(niterations,1);

k = 1;
for i_IMC=1%1:step_size:IMC_steps
    
    disp('--------')
    disp(['IMC step ', num2str(i_IMC)])

    % Load level set
    [a,ed,vp,line_ex,line_ey,line_seg,sep_lines,ngrains,IMC_area,...
    tvec(k),newex,newey,newcoord] = load_level_set(i_IMC,edof_1D); 

    % Store IMC vol/area
    IMC_vol_per_area(k) = IMC_area/L;

    % Plot level set contour 
    f1 = figure(1);
    clf;
    plot_mesh(ex,ey,'black')    
    plot_mesh(newex,newey,'green')
    hold on
    plot_speed_gp(newex,newey,vp)
    for g = 1:ngrains
        g_cols = [2*(g-1) + 1,2*(g-1) + 2];
        plot_interface(line_ex(1:line_seg(g),g_cols),...
                       line_ey(1:line_seg(g),g_cols),my_colors(g,:),g)
    end

    % Axis and title
    axis(axisbc*mm_to_um)
    axis equal
    box on
    title(['Time: ', num2str(tvec(k)), ' (h)'])
    legend

    % Update counter
    k = k + 1;
end

% Save figure
% save_plot(f1,'ls_prof.svg')

% Save IMC volume array
save('mat_eval_files/IMC_volume_MC1.mat','IMC_vol_per_area','tvec')


%% Load functions

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

function [a,ed,vp,line_ex,line_ey,line_seg,sep_lines,ngrains,...
    IMC_area,time,newex,newey,newcoord] ...
    = load_level_set(i_IMC,edof)
% --- Function for loading quantities from a level set step ---

% Filename
s             = what('../single_study/mat_files');
file_location = s.path;
filename      = ['level_set_',num2str(i_IMC),'.mat'];
filename      = fullfile(file_location,filename);

% Load arrays
load(filename,'a')
load(filename,'vpx')
load(filename,'vpy')
load(filename,'line_ex')
load(filename,'line_ey')
load(filename,'line_seg')
load(filename,'sep_lines')
load(filename,'IMC_area')
load(filename,'time')
load(filename,'newex')
load(filename,'newey')
load(filename,'newcoord')


% Rearrange arrays
line_seg   = double(line_seg);
sep_lines  = double(sep_lines);
newex      = newex';
newey      = newey';
newcoord   = newcoord';

% Extract other properties
ngrains = size(a,2);
nelm    = size(edof,1);
ed      = zeros(nelm,4,ngrains);
for g=1:ngrains
    ed(:,:,g) = extract(edof,a(:,g));
end


% Construct vp field from vpx and vpy (3d array output from Fortran not
% possible)
vp        = zeros(2,4,nelm);
vp(1,:,:) = vpx;
vp(2,:,:) = vpy;

% Convert vp to cell struct 
vp1 = repmat({zeros(4,2)},nelm,1);
for ie=1:nelm
    vp1{ie} = vp(:,:,ie)';
end
vp = vp1;

end


function [newex,newey,newcoord] = load_stress(i_IMC)
s = what('../single_study/mat_files');
file_location = s.path;
filename = ['stress_',num2str(i_IMC),'.mat'];
filename = fullfile(file_location,filename);
load(filename,'newex')
load(filename,'newey')
load(filename,'newcoord')

% Rearrange arrays
newcoord = newcoord';
newex    = newex';
newey    = newey'; 

end


function [line_ex,line_ey] = load_unsorted_lines(iload)
s = what('../single_study/mat_files');
file_location = s.path;
filename = ['level_setunsorted_lines_',num2str(iload),'.mat'];
filename = fullfile(file_location,filename);
load(filename,'line_ex')
load(filename,'line_ey')
line_ex = line_ex';
line_ey = line_ey';
end

function vnod = load_vnod_from_diffusion(iload,opt)
% Input data and mesh location
s = what('../single_study/mat_files');
file_location = s.path;
% Input data and mesh location
if (opt==1)
    % IMC phase    
    filename = [file_location, '/diffusion_imc_',num2str(iload),'.mat'];
elseif(opt==2)
    % Sn phase    
    filename = [file_location, '/diffusion_sn_',num2str(iload),'.mat'];
end

% Load data from file
% load(filename,'a')
% load(filename,'r')
% load(filename,'ed')
% load(filename,'j_mean')
load(filename,'vnod')
% ed     = ed';
% j_mean = j_mean';
end

%% Plot functions

function plot_mesh(ex,ey,color)
nm_to_um  = 1e3;

% Plot mesh
patch(ex'*nm_to_um,ey'*nm_to_um,'white','EdgeColor',color,...
    'FaceColor','none','HandleVisibility','off')
end

function plot_interface(line_ex,line_ey,color,g)
mm_to_um  = 1e3;

% % Plot interface
% plot(line_ex'*mm_to_um,line_ey'*mm_to_um,'Color',color,...
%     'LineWidth',2,'Marker','*','DisplayName',['Grain ', num2str(g)])

% Plot interface
plot(line_ex'*mm_to_um,line_ey'*mm_to_um,'Color',...
color,'LineWidth',2,'Marker','o','Markersize',4,'HandleVisibility',...
'off','MarkerFaceColor',color)

plot(line_ex(end,:)*mm_to_um,line_ey(end,:)*mm_to_um,'Color',...
color,'LineWidth',2,'Marker','o','Markersize',4,...
'DisplayName',['GB ' num2str(g)])

end

function plot_speed_ie(ex,ey,vpm)
mm_to_um  = 1e3;
elnods    = size(ex,2);
X         = sum(ex*mm_to_um,2)/elnods;
Y         = sum(ey*mm_to_um,2)/elnods;

% Plot vector v with origin at position (x,y)
sf = 1e4;
quiver(X,Y,vpm(:,1)*sf,vpm(:,2)*sf,'AutoScale','off')

end

function plot_speed_gp(ex,ey,vp)
% Function for plotting velocity in each gauss point

% Change unit to mm
mm_to_um = 1e3;
ex = ex*mm_to_um;
ey = ey*mm_to_um;

% Nelm and ngp
nelm = size(ex,1);
ngp  = size(vp{1},1);

% Allocate quiver plot quantities
X = zeros(nelm*ngp,1);
Y = zeros(nelm*ngp,1);
U = zeros(nelm*ngp,1);
V = zeros(nelm*ngp,1);

% Gauss points
GP   = 0.577350269189626;
xsi4 = [-1  1 -1  1]*GP;
eta4 = [-1 -1  1   1]*GP;

NR4      = zeros(4,4);
NR4(:,1) = (1-xsi4).*(1-eta4);
NR4(:,2) = (1+xsi4).*(1-eta4);
NR4(:,3) = (1+xsi4).*(1+eta4);
NR4(:,4) = (1-xsi4).*(1+eta4);
NR4 = NR4/4;

k = 0; % counter
for ie=1:nelm

    % Element coordinates
    ec = [ex(ie,:)', ey(ie,:)'];
    vp_ie = vp{ie};

    for gp=1:ngp

        % Update counter
        k = k + 1;

        % Global coordinates of gauss point         
        X(k) = NR4(gp,:)*ec(:,1);
        Y(k) = NR4(gp,:)*ec(:,2);

        % Velocity vector at gauss point
        U(k) = vp_ie(gp,1);
        V(k) = vp_ie(gp,2);            
        
    end
end

% Plot vector v with origin at position (x,y)
sf = 1e3;
quiver(X,Y,U*sf,V*sf,'AutoScale','off','Color','k','HandleVisibility','off')
end

function plot_phi_3d(ex,ey,ed,line_ex,line_ey)

% Plot level set function in 3D with fill
patch(ex',ey',ed','white')
view(3)

% Plot level set interface as well
hold on
plot(line_ex',line_ey','k*-')
end

function my_colors = get_my_colors()



% --- Function for providing a set of colors ---

my_colors = [0.8500    0.3250    0.0980
             0.9290    0.6940    0.1250
             0.0046    0.1818    0.3510
             0.6350    0.0780    0.1840
             0         0.4470    0.7410
             0.8687    0.1361    0.0760
             0.0844    0.8693    0.2399
             0.3998    0.5797    0.1233
             0.2599    0.5499    0.1839
             0.8001    0.1450    0.2400];
end


function save_plot(fh,filename)
% Save results in plot dir in single_study if not exist
plotfolder= '../single_study/plots';
if ~exist(plotfolder,"dir")
    mkdir(plotfolder)
    disp(['Created new directory: ', plotfolder])
end    
saveas(fh, fullfile(plotfolder,filename));
end
