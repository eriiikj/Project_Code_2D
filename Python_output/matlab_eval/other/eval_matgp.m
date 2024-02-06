close all; clear; format compact; clc;

% Calfem
addpath /home/er7128ja/Nextcloud/Projekt/Matlab_Calfem/Calfem/fem

%% Load quantities
load 'matout.mat'
load('main_consts.mat')
load('imc_consts.mat')
nelm        = double(nel(1));
nnod        = double(nel(2));
nrgp        = 4;
IMC_y_r     = IMC_r(1);
IMC_x_r1    = IMC_r(2);
IMC_x_r2    = IMC_r(3);
IMC_x_r3    = IMC_r(4);
IMC_x_r4    = IMC_r(5);
grain_width = IMC_r(6);
IMC_steps   = double(IMC_steps(1));
enod        = double(enod);
clear IMC_r nel IMC_steps_in

% Steps
load_steps = double(load_steps);
IMC_steps  = double(IMC_steps);
total_load_steps = load_steps*IMC_steps;
load_step_limits = load_steps*[1:IMC_steps];

right_b_xdofs = (right_boundary_nodes-1)*2 + 1;
left_b_xdofs =  (left_boundary_nodes-1)*2 + 1;

% Convert to standard matlab
ex   = ex';
ey   = ey';

enod = [(1:nelm)', enod'];
ex   = ex';
ey   = ey';

edof = zeros(nelm,9); 
edof(:,1) = enod(:,1);
edof(:,[2,4,6,8]) = (enod(:,[2:end])-1).*2+1;
edof(:,[3,5,7,9]) = (enod(:,[2:end])-1).*2+2;


%% Geometry and unit
% Geometry
model_width  = 4*grain_width;
model_height = max(ey(:));
elm_size     = ex(2,1)-ex(1,1);
nelmx        = model_width/elm_size;

unit = 1e-3; % mm

%% For each iteration, load data and plot
total_load_steps = load_steps*IMC_steps;

figure()
title('IMC radiuses')
hold on
hold on
for i=1:total_load_steps
    load(['matout_iter_', num2str(i),'.mat']);
    IMC_radius = double(IMC_radius);
%     title(['IMC Growth, load step: ', num2str(i)])
    
    plotIMCradius(0,grain_width, IMC_radius   , IMC_x_r1  , IMC_y_r)
    plotIMCradius(grain_width  , grain_width*2, IMC_radius, IMC_x_r2, IMC_y_r)
    plotIMCradius(grain_width*2, grain_width*3, IMC_radius, IMC_x_r3, IMC_y_r)
    plotIMCradius(grain_width*3, grain_width*4, IMC_radius, IMC_x_r4, IMC_y_r)
    axis([0 model_width 1e-3 1.1e-3])
%     eldraw2(ex',ey')
%     for in=1:nnod
%         text(coord(1,in),coord(2,in),num2str(in),'Color','k')
%     end
%     hold on
%     plot_mat_label(nelm,nrgp,gp_mat_label,gpx,gpy)
end

% %% Functions
% function coordout = getcoord(ex,ey,enod)
%   nnod = max(enod(:));
%   coord=zeros(2,nnod);
%   nelm = size(enod,2);
%   ie = 1;
%   for ie=1:nelm
%       coord(1,enod(:,ie)) = ex(:,ie);
%       coord(2,enod(:,ie)) = ey(:,ie);
%   end
%   coordout = coord;
% end
% 
function plotIMCradius(startval,endval, IMC_radius, IMC_x_r, IMC_y_r)
  xvec  = linspace(startval,endval,100);
  xdist = xvec - IMC_x_r;
  ydist = sqrt(IMC_radius^2-xdist.^2);
  yvec  = IMC_y_r-ydist;
  plot(xvec,yvec, 'k','linewidth',1.5)
end
%   
% function plot_mat_label(nelm,nrgp,gp_mat_label,gpx,gpy)
% for ie=1:nelm
%     for j = 1:nrgp
%         if gp_mat_label(j,ie)==1
%             plot(gpx(j,ie),gpy(j,ie),'ko')
%         elseif gp_mat_label(j,ie)==2
%             plot(gpx(j,ie),gpy(j,ie),'ro')
%         elseif gp_mat_label(j,ie)==3
%             plot(gpx(j,ie),gpy(j,ie),'bo')
%         elseif gp_mat_label(j,ie)==4
%             plot(gpx(j,ie),gpy(j,ie),'go')
%         end
%     end
% end
% h = zeros(4, 1);
% h(1) = plot(NaN,NaN,'og');
% h(2) = plot(NaN,NaN,'ob');
% h(3) = plot(NaN,NaN,'or');
% h(4) = plot(NaN,NaN,'ok');
% legend(h, 'Sn','new IMC','old IMC','Cu','AutoUpdate','off');
% end