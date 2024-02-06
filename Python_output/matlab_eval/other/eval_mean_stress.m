clc; close all; clear; format compact; format short

% Path to Calfem
addpath /home/er7128ja/Nextcloud/Projekt/Matlab_Calfem/Calfem/fem

cwd = pwd;

%% Single study
single_path = '/home/er7128ja/Nextcloud/Projekt/Project_Code/Python_output/single_study/';
path_input = single_path;
% single_path = '/home/er7128ja/Nextcloud/Projekt/Project_Code/Python_output/param_study/param_study_mesh_0042';

% --- Display single study ---
% load_and_disp(single_path)

%% Param study

% % Param dirs
% param_path_dir  = '/home/er7128ja/Nextcloud/Projekt/Project_Code/Python_output/param_study/';
% cd(param_path_dir)
% param_dirs     = dir;
% param_dirs     = param_dirs(3:end); % Don't include . and ..
% num_param_dirs = length(param_dirs);
% cd(cwd)
% 
% % --- Display param study ---
% for i=1:num_param_dirs
%     param_dir  = param_dirs(i);
%     param_path = [param_path_dir, param_dir.name];
%     load_and_disp(param_path)
%     
% end


% function load_and_disp(path_input)
    %% Enter path
    cwd = pwd; % Save path
    cd(path_input)
    cd('mat_files')
    pwd
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
    load_step_limits = [1 load_steps*[1:IMC_steps]];

    right_b_xdofs = (right_boundary_nodes-1)*2 + 1;
    left_b_xdofs =  (left_boundary_nodes-1)*2 + 1;

    % Convert to standard matlab
    ex   = ex';
    ey   = ey';
    nelm = size(ex,1);

    enod = [(1:nelm)', enod'];
    ex   = ex';
    ey   = ey';


    edof = zeros(nelm,9); 
    edof(:,1) = enod(:,1);
    edof(:,[2,4,6,8]) = (enod(:,[2:end])-1).*2+1;
    edof(:,[3,5,7,9]) = (enod(:,[2:end])-1).*2+2;


    %% Geometry and unit
    unit = 1e-3; % mm

    % Geometry
    model_width  = 4*grain_width;
    model_height = max(ey(:));
    elm_size     = ey(3,5)-ey(2,5);
    elm_size_y   = 0.025*unit;
    nelmx        = model_width/elm_size;


    

    %% Plot mean stress at given height

    % mean stress heights
    mean_stress_heights = [2.5]*unit;

    % Plot mean stresses
    n_heights = length(mean_stress_heights);
    vm_mean_stress = zeros(total_load_steps,n_heights);
    figure()
    hold on
    box on
    for j=1:n_heights
        mean_stress_height = mean_stress_heights(j);

        a = ey>mean_stress_height;
        a = sum(a,1);
        mean_stress_elms = find(a>0 & a<4);

        % vm mean stress at every load step
        for i=1:total_load_steps
            load([ 'matout_loadstep_', num2str(i),'.mat'])
            vm = rand(size(vm))*10
            vm_mean_stress_elms   = vm(:,mean_stress_elms);
            vm_mean_stress(i,j)   = sum(vm_mean_stress_elms(:))/...
                                length(vm_mean_stress_elms(:));
        end
        plot(1:total_load_steps,vm_mean_stress(:,j), '--*')
    end

    for i=1:length(load_step_limits)-1
        xline(load_step_limits(i),'--',['IMC step ',num2str(i)],...
              'LabelVerticalAlignment','middle');
    end
    axis([1-0.2 total_load_steps+1 0 13])
    title(['Mean von Mises stress. ',...
           'Volume inc: ', num2str(IMC_vol_inc), '. ',...
           'Mesh size: ', sprintf('%.3f',elm_size*1e3), ' microns.'], ...
            'Interpreter','latex')
    xlabel('Load steps')
    ylabel('Mean von Mises stress [MPa]')

    legend('1.1','1.2','1.5')
    
    %% Leave dir
%     cd(cwd)

% end
