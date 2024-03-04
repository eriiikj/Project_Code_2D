% -----------------------------------
% Extracting and plotting diffusion phase properties  
% from .mat files generated by Fortran program.
% Written by Erik Jacobsson 2023-04-12
% -----------------------------------

clear; clc; close all;

%% Define the grain struct

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

% Import level set mesh properties
[ngrains,IMC_steps,enod,nodel,nelm,nnod,edof] = load_level_set_init();


% Create an array to store multiple grains
grainArr = repmat(grain, 1, ngrains);

%% Diffusion equation
step_size = 1;
for i_IMC = 1%1:step_size:IMC_steps


    % Load level set
    [a,ed,line_ex,line_ey,line_seg,line_coord,line_coordN,ex,ey,coord,...
        tppoints] = load_level_set(i_IMC,edof);

    % Load triangle mesh
    [coordT,enodT,edofT,nelmT,ndofT,dofspernodeT,nodelT,exT,eyT] = ...
    import_triangle_mesh(i_IMC);

    % Load interpolated level set field
    [ls,ls_ed] = import_diffusion_glob(i_IMC,edofT,nodelT);

    figure(13)
    g = 2;
%     plot_field2(ex,ey,ed(:,:,g),'k')
    hold on
    plot_field(exT,eyT,ls_ed(:,:,g),'k')    
    g_cols = [2*(g-1) + 1,2*(g-1) + 2];
    plot_interface(line_ex(1:line_seg(g),g_cols),...
                   line_ey(1:line_seg(g),g_cols),'k',g)
    axis equal
    box on


    % Extract nodes
    elmsTG  = sum((ls_ed(:,:,g)<=0),2)==nodelT;
    nelmTG  = sum(elmsTG)
    enodTG  = enodT(elmsTG,:);   
    nodsTG  = unique(enodTG);
    coordTG = coordT(nodsTG,:);
    nnodTG  = size(coordTG,1)
    

    % New enod
    newGrainNods = (1:nnodTG)';
    for k=1:nnodTG
        mask = enodTG == nodsTG(k);
        enodTG(mask) = newGrainNods(k);
    end
    edofTG = [(1:nelmTG)',enodTG];
    dofsTG = (1:nnodTG)';

    % Ex and Ey
    [exTG,eyTG] = coordxtr(edofTG,coordTG,dofsTG,nodelT);

    figure()
    patch(exTG',eyTG','white')
    axis equal
    box on

end



%% Load functions

function [ngrains,IMC_steps,enod,nodel,nelm,nnod,edof_1D] = ...
    load_level_set_init()
% --- Function for loading initial level set data ---

% Filename
s             = what('../single_study/mat_files');
file_location = s.path;
filename      = 'level_setinit.mat';
filename      = fullfile(file_location,filename);

% Load arrays
load(filename,'enod')
load(filename,'nodel')
load(filename,'ngrains')

% Rearrange arrays
enod             = double(enod)';
nodel            = double(nodel);

% Extract other properties
nelm             = size(enod,1);
nnod             = max(enod(:));
edof_1D          = zeros(nelm,nodel+1);
edof_1D(:,1)     = 1:nelm;
edof_1D(:,2:end) = enod;

% IMC steps (finished level set updates)
filename_prefix = 'level_set_';
files           = dir(fullfile(file_location, [filename_prefix, '*']));
IMC_steps       = numel(files);

end


function [a,ed,line_ex,line_ey,line_seg,line_coord,line_coordN,...
    newex,newey,newcoord,tppoints] = load_level_set(i_IMC,edof)
% --- Function for loading quantities from a level set step ---

% Filename
s             = what('../single_study/mat_files');
file_location = s.path;
filename      = ['level_set_',num2str(i_IMC),'.mat'];
filename      = fullfile(file_location,filename);

% Load arrays
load(filename,'a')
load(filename,'line_ex')
load(filename,'line_ey')
load(filename,'line_seg')
load(filename,'newex')
load(filename,'newey')
load(filename,'newcoord')
load(filename,'line_coord')
load(filename,'line_coordN')
load(filename,'tppoints')

% Rearrange arrays
line_seg    = double(line_seg);
line_coordN = double(line_coordN);
newex       = newex';
newey       = newey';
newcoord    = newcoord';

ngrains = size(a,2);
nelm    = size(edof,1);
ed      = zeros(nelm,4,ngrains);
for g=1:ngrains
    ed(:,:,g) = extract(edof,a(:,g));
end

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


function [coord,enod,edof,nelm,ndof,dofspernode,nodel,ex,ey] = ...
    import_triangle_mesh(i_IMC)


% Input data and mesh location
s = what('../single_study/phase_meshes');
input_mesh_location = s.path;
filename=[input_mesh_location, '/triangle_mesh_',num2str(i_IMC),'.json'];

% Import mesh_struct from json file
fid = fopen(filename); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
mesh_struct = jsondecode(str);
clear fid raw str

% Extract mesh data from mesh_struct
enod           = mesh_struct.enod';
coord          = mesh_struct.coord'*1e-3;
nodel          = mesh_struct.nodel;
nelm           = mesh_struct.nelm;
ndof           = mesh_struct.ndof;
dofspernode    = mesh_struct.dofs_per_node;

% Dofs (one dof per node)
if (dofspernode==1)
    dofs = (1:ndof)';
end

% Edof
edof = [(1:nelm)',enod];

% Ex and Ey
[ex,ey] = coordxtr(edof,coord,dofs,nodel);
end

function [ls,ls_ed] = import_diffusion_glob(i_IMC, edof, nodel)
% Input data and mesh location
s = what('../single_study/mat_files/');
input_mesh_location = s.path;
filename=[input_mesh_location, '/triangle_ls_',num2str(i_IMC),'.mat'];
load(filename,'ls')

ngrains = size(ls,2);
nelm    = size(edof,1);
ls_ed      = zeros(nelm,nodel,ngrains);
for g=1:ngrains
    ls_ed(:,:,g) = extract(edof,ls(:,g));
end

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



function coord = lines_to_coord(line_ex, line_ey)
    % Function for extracting unique vertices in lines
    N = size(line_ex, 1);
    
    % Initialize coord matrix
    coord = zeros(2 * N, 2);
    
    % Assign values to coord matrix
    coord(1:N, 1) = line_ex(:, 1);
    coord(N+1:end, 1) = line_ex(:, 2);
    coord(1:N, 2) = line_ey(:, 1);
    coord(N+1:end, 2) = line_ey(:, 2);
    
    % Find unique vertices
    unique_coord = unique(coord, 'rows');
    
    % Return coord matrix
    coord = unique_coord;
end


%% Plot functions

function plot_mesh(ex,ey,color)
nm_to_um  = 1e3;

% Plot mesh
patch(ex'*nm_to_um,ey'*nm_to_um,'white','EdgeColor',color,...
    'FaceColor','none','HandleVisibility','off')
end

function plot_field(ex,ey,ed,color)
nm_to_um  = 1e3;

% Plot mesh
patch(ex'*nm_to_um,ey'*nm_to_um,ed',ed','EdgeColor',color,...
    'HandleVisibility','off')

cc                = colorbar;
cc.Label.String   = '\bf Value';
cc.Label.FontSize = 10;

end

function plot_field2(ex,ey,ed,color)
nm_to_um  = 1e3;

% Plot mesh
patch(ex'*nm_to_um,ey'*nm_to_um,ed','white','EdgeColor',color,...
    'FaceColor','none','HandleVisibility','off')

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


function plot_interface(line_ex,line_ey,color,g)
mm_to_um  = 1e3;

% % Plot interface
% plot(line_ex'*mm_to_um,line_ey'*mm_to_um,'Color',color,...
%     'LineWidth',2,'Marker','*','DisplayName',['Grain ', num2str(g)])

% Plot interface
plot(line_ex'*mm_to_um,line_ey'*mm_to_um,'Color',...
color,'LineWidth',3,'Marker','o','Markersize',4,'HandleVisibility',...
'off','MarkerFaceColor',color)

plot(line_ex(end,:)*mm_to_um,line_ey(end,:)*mm_to_um,'Color',...
color,'LineWidth',3,'Marker','o','Markersize',4,...
'DisplayName',['GB ' num2str(g)])

end
