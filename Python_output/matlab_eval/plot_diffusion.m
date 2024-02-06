% -----------------------------------
% Extracting and plotting diffusion phase properties  
% from .mat files generated by Fortran program.
% Written by Erik Jacobsson 2023-04-12
% -----------------------------------

clear; clc; close all;

%% Define the grain struct

% Load level set init
[coord,enod,nodel,ex,ey,gpx,gpy,nelm,nnod,edof_1D,axisbc,IMC_steps,...
    ngrains] = load_level_set_init();

% Define the grain struct
grain = struct(...
    'edof', [], ...     
    'enod', [], ...     
    'coord', [], ...    
    'dofs', [], ...     
    'ex', [], ...       
    'ey', [], ...      
    'bcnod', [],...
    'bcval',[],...
    'bcval_idx',[],...    
    'nodel',[],...
    'nelm',[],...
    'ndof',[],...
    'a',[],...
    'r',[],...
    'ed',[],...
    'j_flux',[],...
    'vnod',[],...
    'p',[],...
    'p_ed',[]);

% Create an array to store multiple grains
grainArr = repmat(grain, 1, ngrains);

%% Diffusion equation
step_size = 1;
for i_IMC = 1%1:step_size:IMC_steps

    % Load ngrains from level set
    [ngrains] = load_level_set(i_IMC); 
    
    % Load grain mesh
    for g=1:ngrains
        [grainArr(g).edof,grainArr(g).enod,grainArr(g).coord,...
         grainArr(g).dofs,grainArr(g).ex,grainArr(g).ey,...
         grainArr(g).bcnod,grainArr(g).bcval,grainArr(g).bcval_idx,...
         grainArr(g).nodel,grainArr(g).nelm,grainArr(g).ndof]...
         = import_grain_mesh(i_IMC,g);
    
        % Load grain results
        [grainArr(g).a,grainArr(g).r,grainArr(g).ed,...
         grainArr(g).jint, grainArr(g).j_flux,grainArr(g).p,...
         grainArr(g).p_ed] = load_diffusion(i_IMC,g,grainArr(g).edof);
    end


    % Plot grain mesh with concentration
    f1 = figure(1);
    cla;
    for g=1:ngrains
%         plot_mesh(grainArr(g).ex,grainArr(g).ey)
        plot_2D_conc(grainArr(g).ex,grainArr(g).ey,grainArr(g).ed,i_IMC)
        hold on
        plot_j(grainArr(g).ex,grainArr(g).ey,grainArr(g).j_flux)
        axis equal  
        box on
    end

% %     Plot grain mesh with hydrostatic pressure
%     f2 = figure(2);
%     cla;
%     for g=4
%         plot_mesh(grainArr(g).ex,grainArr(g).ey)
%         plot_2D_pressure(grainArr(g).ex,grainArr(g).ey,...
%         grainArr(g).p_ed,i_IMC)
%         hold on
%         axis equal
%         axis(axisbc*1e3)    
%         box on
%         clim([-10 10])
%     end
end

% Save figure
% save_plot(f1,'conc_prof.svg')



% % --- Compute jint2 ---
% g = 4;
% nbnods  = length(grainArr(g).bcnod(:,1));
% jint2   = zeros(nbnods,1);
% counter = 1;
% while (counter<nbnods)
%     iseg_start  = counter;
%     bcval_c     = grainArr(g).bcval_idx(counter);
%     nnods_bcseg = 0;   
%     stop        = 0;
%     while (stop==0 && (counter + nnods_bcseg<=nbnods))        
%         bcval_cc    = grainArr(g).bcval_idx(counter + nnods_bcseg);
%         if (abs(bcval_cc-bcval_c)<1e-15) 
%             nnods_bcseg = nnods_bcseg + 1;
%         else
%             stop = 1;
%         end
%     end
%     iseg_end    = iseg_start + nnods_bcseg - 1;
%     nseg_sep    = nnods_bcseg - 1;
%     if (nnods_bcseg>1)
%    
%         % Compute mass matrix
%         enod_m = zeros(2,nseg_sep);
%         for ie=1:nseg_sep
%           enod_m(:,ie) = [ie,ie+1];
%         end
%         M = zeros(nnods_bcseg,nnods_bcseg);
%         
%         % Assemble mass matrix
%         for ie=1:nseg_sep
%         
%           % Nods
%           n1 = grainArr(g).bcnod(counter+ie-1,1);
%           n2 = grainArr(g).bcnod(counter+ie,1); 
%         
%           % Coordinates and length of element
%           x1 = grainArr(g).coord(n1,1);
%           y1 = grainArr(g).coord(n1,2);
%           x2 = grainArr(g).coord(n2,1);
%           y2 = grainArr(g).coord(n2,2);
%           Le = norm([(x2-x1), (y2-y1)]);
%         
%           % Mass matrix
%           Me      = zeros(2,2);
%           Me(1,:) = [2, 1];
%           Me(2,:) = [1, 2];
%           Me      = Le/6d0*Me;
%           M(enod_m(:,ie),enod_m(:,ie)) = M(enod_m(:,ie),enod_m(:,ie)) + Me;
%         end    
%     
%         % Solve M*vnod = -r
%         b = -grainArr(g).r(grainArr(g).bcnod(iseg_start:iseg_end,1));
%         jint2(iseg_start:iseg_end) = M\b;
% 
%     end
% 
%     %  Update counter
%     counter = counter + nnods_bcseg;
% 
% end



%% Load and extract functions

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


function [edof,enod,coord,dofs,ex,ey,bcnod,bcval,bcval_idx,nodel, ...
    nelm,ndof] = import_grain_mesh(i_IMC,g)


% Input data and mesh location
s = what('../single_study/phase_meshes');
input_mesh_location = s.path;
filename=[input_mesh_location, '/phase_mesh_',num2str(i_IMC),'_',...
          num2str(g),'.json'];

% Extracting mesh from python created json file
fid = fopen(filename); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
mesh_struct = jsondecode(str);
clear fid raw str
edof           = zeros(mesh_struct.nelm, mesh_struct.dofel+1);
edof(:,1)      = 1:mesh_struct.nelm;
edof(:,2:end)  = mesh_struct.edof;
enod           = mesh_struct.enod';
bcnod          = mesh_struct.bcnod;
bcval          = mesh_struct.bcval;
bcval_idx      = mesh_struct.bcval_idx;
coord          = mesh_struct.coord'*1e-3;
nodel          = mesh_struct.nodel;
nelm           = mesh_struct.nelm;
ndof           = mesh_struct.ndof;
dofspernode    = mesh_struct.dofs_per_node;

% Dofs (one dof per node)
if (dofspernode==1)
    dofs = (1:ndof)';
end

% Ex and Ey
[ex,ey] = coordxtr(edof,coord,dofs,nodel);
end


function [a,r,ed,jint,j_flux,p,p_ed] = load_diffusion(i_IMC,g,edof)

% Input data and mesh location
s = what('../single_study/mat_files');
file_location = s.path;
% Input data and mesh location
filename = [file_location, '/diffusion_',num2str(i_IMC),'_',...
    num2str(g),'.mat'];

% Load data from file
load(filename,'a')
load(filename,'r')
load(filename,'ed')
load(filename,'jint')
load(filename,'j_flux')
load(filename,'p')

% Reshape arrays
ed     = ed';
j_flux = j_flux';
if (length(p)>1)
    p_ed = extract(edof,p);
else
    p_ed = 0;
end

end


function [ngrains] = load_level_set(i_IMC)
% --- Function for loading quantities from a level set step ---

% Filename
s             = what('../single_study/mat_files');
file_location = s.path;
filename      = ['level_set_',num2str(i_IMC),'.mat'];
filename      = fullfile(file_location,filename);

% Load arrays
load(filename,'a')

% Extract other properties
ngrains = size(a,2);

end


%% Plot functions

function plot_mesh(ex,ey)
nm_to_um  = 1e3;

% Plot mesh
patch(ex'*nm_to_um,ey'*nm_to_um,'white')
end

function plot_2D_conc(ex,ey,ed,i_IMC)
nm_to_um = 1e3;

p = patch(ex'*nm_to_um,ey'*nm_to_um,ed');
% p.EdgeColor = 'interp';

cc                = colorbar;
cc.Label.String   = '\bf Diffusion potential (mol/mm^3)';
cc.Label.FontSize = 10;
xlabel('\bf Micrometer')
ylabel('\bf Micrometer')
title(['IMC step: ', num2str(i_IMC)])
end

function plot_2D_pressure(ex,ey,ed,i_IMC)
nm_to_um = 1e3;

p = patch(ex'*nm_to_um,ey'*nm_to_um,ed');
% p.EdgeColor = 'interp';

cc                = colorbar;
cc.Label.String   = '\bf Hydrostatic pressure (MPa)';
cc.Label.FontSize = 10;
xlabel('\bf Micrometer')
ylabel('\bf Micrometer')
title(['IMC step: ', num2str(i_IMC)])
end



function plot_j(ex,ey,j)
nm_to_um = 1e3;
elnods = size(ex,2);
X = sum(ex*nm_to_um,2)/elnods;
Y = sum(ey*nm_to_um,2)/elnods;

% Plot vector field
sf = 1e8;
quiver(X,Y,j(:,1)*sf,j(:,2)*sf,'AutoScale','off','Color','k')
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
