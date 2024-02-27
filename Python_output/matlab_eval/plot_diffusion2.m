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
    'nnod',[],...
    'a',[],...
    'r',[],...
    'ed',[],...
    'j_flux',[],...
    'vnod',[],...
    'p',[],...
    'p_ed',[]);

% Import level set mesh properties
[ngrains,IMC_steps,~,~,~,~,edof] = load_level_set_init();


% Create an array to store multiple grains
grainArr = repmat(grain, 1, ngrains);

%% Diffusion equation
step_size = 1;
for i_IMC = 1%1:step_size:IMC_steps

    % Load level set
    [a,ed,line_ex,line_ey,line_seg,line_coord,line_coordN,ex,ey,coord,...
        tppoints] = load_level_set(i_IMC,edof);

    % Load grain mesh
    for g=1:ngrains
        [grainArr(g).edof,grainArr(g).enod,grainArr(g).coord,...
         grainArr(g).dofs,grainArr(g).ex,grainArr(g).ey,...
         grainArr(g).bcnod,grainArr(g).bcval,...
         grainArr(g).nodel,grainArr(g).nelm,grainArr(g).nnod]...
         = import_grain_mesh(i_IMC,g);
    
        % Load grain results
        [grainArr(g).a,grainArr(g).r,grainArr(g).ed,...
         grainArr(g).jint, grainArr(g).j_flux,grainArr(g).p,...
         grainArr(g).p_ed] = load_diffusion(i_IMC,g,grainArr(g).edof);
    end


    % Plot grain mesh with concentration
    f1 = figure(1);
    cla;
    for g=1%1:ngrains
%         plot_mesh(grainArr(g).ex,grainArr(g).ey,'k')
        plot_2D_conc(grainArr(g).ex,grainArr(g).ey,grainArr(g).ed,i_IMC)
        hold on
        plot_j(grainArr(g).ex,grainArr(g).ey,grainArr(g).j_flux)
        axis equal  
        box on
    end

% %     Plot grain mesh with hydrostatic pressure
%     f2 = figure(2);
%     cla;
%     for g=1:ngrains
%         plot_mesh(grainArr(g).ex,grainArr(g).ey)
%         plot_2D_pressure(grainArr(g).ex,grainArr(g).ey,...
%         grainArr(g).p_ed,i_IMC)
%         hold on
%         axis equal 
%         box on
%         clim([-10 10])
%     end






%     % Load level set
%     [a,ed,line_ex,line_ey,line_seg,line_coord,line_coordN,ex,ey,coord,...
%         tppoints] = load_level_set(i_IMC,edof);
% 
%     % Load triangle mesh
%     [coordT,enodT,edofT,nelmT,ndofT,dofspernodeT,nodelT,exT,eyT] = ...
%     import_triangle_mesh(i_IMC);
% 
%     % Load interpolated level set field
%     [ls,ls_ed] = import_diffusion_glob(i_IMC,edofT,nodelT);


%     figure(13)
%     g = 1;
%     plot_field2(ex,ey,ed(:,:,g),'k')
%     hold on
%     plot_field(exT,eyT,ls_ed(:,:,g),'k')    
%     g_cols = [2*(g-1) + 1,2*(g-1) + 2];
%     plot_interface(line_ex(1:line_seg(g),g_cols),...
%                    line_ey(1:line_seg(g),g_cols),'k',g)
%     axis equal
%     box on
% 
% 
%     % Extract nodes
%     elmsTG  = sum((ls_ed(:,:,g)<=0),2)==nodelT;
%     nelmTG  = sum(elmsTG);
%     enodTG  = enodT(elmsTG,:);   
%     nodsTG  = unique(enodTG);
%     coordTG = coordT(nodsTG,:);
%     nnodTG  = size(coordTG,1);
% %     exTG    = exT(elmsTG,:);
% %     eyTG    = eyT(elmsTG,:);
%     
% 
%     % New enod
%     newGrainNods = (1:nnodTG)';
%     for k=1:nnodTG
%         mask = enodTG == nodsTG(k);
%         enodTG(mask) = newGrainNods(k);
%     end
%     edofTG = [(1:nelmTG)',enodTG];
%     dofsTG = (1:nnodTG)';
% 
%     % Ex and Ey
%     [exTG,eyTG] = coordxtr(edofTG,coordTG,dofsTG,nodelT);
% 
%     figure()
%     patch(exTG',eyTG','white')
%     axis equal
%     box on



end


g_cols  = [2*(g-1) + 1,2*(g-1) + 2];
nseg = line_seg(g);
line_ex = line_ex(1:nseg,g_cols);
line_ey = line_ey(1:nseg,g_cols);


% --- Compute jint2 ---
nbnods  = length(grainArr(g).bcnod(:,1));
jint2   = zeros(nbnods,1);
bcvals  = grainArr(g).bcval;
bcnods  = grainArr(g).bcnod(:,1);
nbcnods = length(bcnods);
bcval_u = unique(bcvals);
% Loop through each bcval (phase interface). 
% 1) Find coordinates along the bcval interface
% 2) Find the lines connecting the coordinates
bool_seg = false(nseg, 1);
for j = 1:length(bcval_u)
    bool_seg(:)  = false;
    bcval        = bcval_u(j);
    bcvals_b     = (bcvals == bcval);
    nnods_bcseg  = sum(bcvals_b);
    bcval_coords = grainArr(g).coord(bcnods(bcvals_b), :);
    enodL        = zeros(1000,2);
    cc = 0;
    for iseg = 1:nseg
        x1 = line_ex(iseg, 1); x2 = line_ex(iseg, 2);
        y1 = line_ey(iseg, 1); y2 = line_ey(iseg, 2);
        
        P1incoords = any((abs(bcval_coords(:,1) - x1) < 1e-8) &...
                          abs(bcval_coords(:,2) - y1) < 1e-8);
        P2incoords = any((abs(bcval_coords(:,1) - x2) < 1e-8) &...
                          abs(bcval_coords(:,2) - y2) < 1e-8);
        if (P1incoords && P2incoords)
            bool_seg(iseg) = true;
            nod1 = find(abs(grainArr(g).coord(:,1)-x1)<1e-8 & abs(grainArr(g).coord(:,2)-y1)<1e-8);
            nod2 = find(abs(grainArr(g).coord(:,1)-x2)<1e-8 & abs(grainArr(g).coord(:,2)-y2)<1e-8);
            cc =  cc + 1;
            enodL(cc,:) = [nod1 nod2];
            
        end
    end 
    enodL      = enodL(1:cc,:);
    % New enod
    oldNods = unique(enodL(:));
    newNods = (1:length(oldNods))';
    for k=1:length(oldNods)
        mask = enodL == oldNods(k);
        enodL(mask) = newNods(k);
    end
    enod_m     = enodL';
    lines      = find(bool_seg);
    nseg_bcval = length(lines);
    M = zeros(length(oldNods),length(oldNods));
    for k=1:nseg_bcval
        iseg = lines(k);
        x1 = line_ex(iseg, 1); x2 = line_ex(iseg, 2);
        y1 = line_ey(iseg, 1); y2 = line_ey(iseg, 2);
        Le = norm([(x2-x1), (y2-y1)]);
        % Mass matrix
        Me      = zeros(2,2);
        Me(1,:) = [2, 1];
        Me(2,:) = [1, 2];
        Me      = Le/6d0*Me;
        M(enod_m(:,k),enod_m(:,k)) = M(enod_m(:,k),enod_m(:,k)) + Me;
    end

    % Solve M*vnod = -r
    b = -grainArr(g).r(bcnods(bcvals_b));
    jint2(bcvals_b) = M\b;
end

% while (counter<nbnods)
%     iseg_start  = counter;
%     bcval_c     = grainArr(g).bcval(counter);
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


function [edof,enod,coord,dofs,ex,ey,bcnod,bcval,nodel, ...
    nelm,nnod] = import_grain_mesh(i_IMC,g)


% Input data and mesh location
s = what('../single_study/phase_meshes');
input_mesh_location = s.path;
filename=[input_mesh_location, '/phase_mesh_',num2str(i_IMC),'_g',...
          num2str(g),'.json'];

% Extracting mesh from python created json file
fid = fopen(filename); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
mesh_struct = jsondecode(str);
clear fid raw str
edof           = zeros(mesh_struct.nelm, mesh_struct.nodel+1);
edof(:,1)      = 1:mesh_struct.nelm;
edof(:,2:end)  = mesh_struct.enod';
enod           = mesh_struct.enod';
bcnod          = mesh_struct.bcnod;
bcval          = mesh_struct.bcval;
% bcval_idx      = mesh_struct.bcval_idx;
coord          = mesh_struct.coord'*1e-3;
nodel          = mesh_struct.nodel;
nelm           = mesh_struct.nelm;
nnod           = mesh_struct.nnod;

% Dofs (one dof per node)
dofs = (1:nnod)';

% Ex and Ey
[ex,ey] = coordxtr(edof,coord,dofs,nodel);
end


function [coord,enod,edof,nelm,ndof,dofspernode,nodel,ex,ey] = ...
    import_triangle_mesh(i_IMC)


% Input data and mesh location
s = what('../single_study/phase_meshes');
input_mesh_location = s.path;
filename=[input_mesh_location, '/phase_mesh_',num2str(i_IMC),'.json'];

% Import mesh_struct from json file
fid = fopen(filename); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
mesh_struct = jsondecode(str);
clear fid raw str

% Extract mesh data from mesh_struct
enod           = mesh_struct.enod' + 1;
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
filename=[input_mesh_location, '/diffusion_',num2str(i_IMC),'.mat'];
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