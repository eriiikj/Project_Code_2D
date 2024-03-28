% -----------------------------------
% Extracting and plotting stresses from
% .mat files generated by Fortran program.
% Written by Erik Jacobsson 2023-02-16
% -----------------------------------

% clear; clc; close all;
close all; clear;


%% Import data

% Load level set init
[coord,enod,nodel,ex,ey,gpx,gpy,nelm,nnod,edof_1D,axisbc,IMC_steps,...
    ngrains] = load_level_set_init();

% Length
L = axisbc(2)-axisbc(1);

% Load initial morphology
i_IMC=1;
[a,ed,vp,line_ex_1,line_ey_1,line_seg_1,sep_lines_1,ngrains,IMC_area,...
t] = load_level_set(i_IMC,edof_1D);


%% Load loop
IMC_steps   = 190;
step_size   = 5;
niterations = floor((IMC_steps - 1) / step_size);

% Storing all times
tvec             = zeros(niterations,1);
IMC_vol_per_area = zeros(niterations,1);

k = 1;
for i_IMC=35%1:step_size:IMC_steps

    disp('--------')
    disp(['IMC step ', num2str(i_IMC)])  
    
    % Load time from level set
    [~,~,~,~,~,~,~,~,~,~,time] = load_level_set(i_IMC,edof_1D);

    % Load level set
    [a,ed,vp,line_ex,line_ey,line_seg,sep_lines,ngrains,IMC_area,...
    hphi_ed,tvec(k)] = load_level_set(i_IMC,edof_1D);

    % Store IMC vol/area
    IMC_vol_per_area(k) = IMC_area/L;

    % Load stress state
    [vm,vm_Cu,vm_IMC,vm_Sn,biax,biax_Cu,biax_IMC,biax_Sn,biax_Sn2,p_ed,...
     newcoord,newex,newey] = load_stress(i_IMC,edof_1D);    

    axisbc = [min(newcoord(:,1)),max(newcoord(:,1)),min(newcoord(:,2)),max(newcoord(:,2))];
%     axisbc(2) = round((axisbc(2)*1e3),2)*1e-3;
    disp(['Time: ', num2str(tvec(k))])

%     %  --- von Mises stress ---
%     f1 = figure(1);
%     cla;    
%     plot_2D_stress(newex,newey,vm,axisbc,...
%         'von Mises stress (MPa)');
%     hold on
% %     for g = 2:ngrains-5
% %         g_cols = [2*(g-1) + 1,2*(g-1) + 2];
% %         plot_interface(line_ex_1(1:line_seg(g),g_cols),...
% %                        line_ey_1(1:line_seg(g),g_cols),'w',':')      
% %     end
%     for g = 1:ngrains
%         g_cols = [2*(g-1) + 1,2*(g-1) + 2];
%         plot_interface(line_ex(1:line_seg(g),g_cols),...
%                        line_ey(1:line_seg(g),g_cols),'w','-')        
%     end
%     % Save
%     plotfolder = '../single_study/plots';
%     filename   = 'vm2.eps';    
%     exportgraphics(gca,fullfile(plotfolder,filename),'Resolution',300)
 
    %  --- von Mises stress in Sn ---
    f2 = figure(2);
    cla;    
    plot_2D_stress(newex,newey,vm_Sn,axisbc,...
        'von Mises stress (MPa)');
    hold on
%     for g = 2:ngrains-5
%         g_cols = [2*(g-1) + 1,2*(g-1) + 2];
%         plot_interface(line_ex_1(1:line_seg(g),g_cols),...
%                        line_ey_1(1:line_seg(g),g_cols),'w',':')      
%     end
    for g = 1:ngrains
        g_cols = [2*(g-1) + 1,2*(g-1) + 2];
        plot_interface(line_ex(1:line_seg(g),g_cols),...
                       line_ey(1:line_seg(g),g_cols),'w','-')        
    end
    clim([0 14.5])
    % Save
    plotfolder = '../single_study/plots';
    filename   = 'vmSn1.png';    
    exportgraphics(gca,fullfile(plotfolder,filename),'Resolution',300);
%     saveas(gca,fullfile(plotfolder,filename));


% %  --- Biaxial stress in Sn 2 ---
%     f4 = figure(4);
%     cla;    
%     plot_2D_stress(newex,newey,biax_Sn2,axisbc,...
%         'Biaxial stress (MPa)');
%     hold on
%     for g = 2:ngrains-5
%         g_cols = [2*(g-1) + 1,2*(g-1) + 2];
%         plot_interface(line_ex_1(1:line_seg(g),g_cols),...
%                        line_ey_1(1:line_seg(g),g_cols),'w',':')      
%     end
%     for g = 1:ngrains
%         g_cols = [2*(g-1) + 1,2*(g-1) + 2];
%         plot_interface(line_ex(1:line_seg(g),g_cols),...
%                        line_ey(1:line_seg(g),g_cols),'w','-')        
%     end
%     clim([-25 15])


    % --- Volume expansion ---
%     figure()
% %     hp = 0.45*max(max(hphi_ed(:,:,2),hphi_ed(:,:,3)),hphi_ed(:,:,4));
% %     hp = 0.45*max(hphi_ed(:,:,2),hphi_ed(:,:,3));
%     hp = 0.45*max(max(max(hphi_ed(:,:,2),hphi_ed(:,:,3)),hphi_ed(:,:,4)),hphi_ed(:,:,5));
%     plot_2D_stress(newex,newey,hp,axisbc,...
%         'Volume expansion');

    % Update counter
    k = k + 1;
end

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

function [vm,vm_Cu,vm_IMC,vm_Sn,biax,biax_Cu,biax_IMC,biax_Sn,biax_Sn2,...
    p_ed,newcoord,newex,newey] = load_stress(i_IMC,edof_1D)
s = what('../single_study/mat_files');
file_location = s.path;
filename = ['stress_',num2str(i_IMC),'.mat'];
filename = fullfile(file_location,filename);
load(filename,'vm_nodplot')
load(filename,'vm_Cu_nodplot')
load(filename,'vm_IMC_nodplot')
load(filename,'vm_Sn_nodplot')
load(filename,'biax_nodplot')
load(filename,'biax_Cu_nodplot')
load(filename,'biax_IMC_nodplot')
load(filename,'biax_Sn_nodplot')
load(filename,'biax_Sn2_nodplot')
load(filename,'newex')
load(filename,'newey')
load(filename,'newcoord')

% Rearrange arrays
newcoord = newcoord';
newex    = newex';
newey    = newey'; 

% Hydrostatic pressure
load(filename,'p_nod')

% Extract
vm       = extract(edof_1D,vm_nodplot);
vm_Cu    = extract(edof_1D,vm_Cu_nodplot);
vm_IMC   = extract(edof_1D,vm_IMC_nodplot);
vm_Sn    = extract(edof_1D,vm_Sn_nodplot);
biax     = extract(edof_1D,biax_nodplot);
biax_Cu  = extract(edof_1D,biax_Cu_nodplot);
biax_IMC = extract(edof_1D,biax_IMC_nodplot);
biax_Sn  = extract(edof_1D,biax_Sn_nodplot);
biax_Sn2 = extract(edof_1D,biax_Sn2_nodplot);
p_ed     = extract(edof_1D,p_nod);
end

function [cc] = plot_2D_stress(ex,ey,ed,axisbc,stress_str)
mm_to_um = 1e3;

p           = patch(ex'*mm_to_um,ey'*mm_to_um,ed');
p.EdgeColor = 'interp';

% Colormap
colormap("parula");
% cc                      = colorbar;
% cctitle                 = stress_str;
% cc.Label.String         = ['\bf ', cctitle];
% cc.Label.Interpreter    = 'latex';
% cc.Ruler.Exponent       = 0;
% cc.TickLabelInterpreter = 'latex';
% cc.FontSize             = 16;
% cc.Label.FontSize       = 16;
% set(cc,'YTick',[0,2,4,6,8,10,12,14])

% Labels and axis
xlabel('\textbf{Micrometer}','Interpreter','Latex')
ylabel('\textbf{Micrometer}','Interpreter','Latex')
% title([' Time: ', num2str(i_IMC), ' (h). Material: ', mat_str, '.'])
axis equal
axis(axisbc*mm_to_um)    
box on
set(gca,'Yticklabel',[])
% ytickformat('%.2f')
% xticks([axisbc(1)*mm_to_um,1,2,3,4,axisbc(2)*mm_to_um])
xticks([axisbc(1)*mm_to_um,0.5,axisbc(2)*mm_to_um])
xtickformat('%.2f')
ax = gca;
ax.FontSize = 16;
set(gca, 'TickLabelInterpreter', 'latex');
box on
ax.Position = [0.1120 0.1482 0.6676 0.7768];
% ax.Position = [0.1120 0.0300 0.6676 0.8150];
% ax.Position = [0.1123 0.1512 0.6695 0.7738];
set(gca, 'TickDir', 'both', 'Layer', 'top');
end



function [a,ed,vp,line_ex,line_ey,line_seg,sep_lines,ngrains,...
    IMC_area,hphi_ed,time] = load_level_set(i_IMC,edof)
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
load(filename,'hphi')
load(filename,'time')

% Rearrange arrays
line_seg  = double(line_seg);
sep_lines = double(sep_lines);


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

hphi_ed  = zeros(nelm,4,ngrains);
for g=1:ngrains
    hphi_ed(:,:,g) = extract(edof,hphi(:,g));
end
end


function plot_interface(line_ex,line_ey,color,marker)
mm_to_um  = 1e3;

% Plot interface
plot(line_ex'*mm_to_um,line_ey'*mm_to_um,'Color',color,...
    'LineWidth',2,'LineStyle',marker)
end
