% -----------------------------------
% Extracting and plotting times from .mat files generated
% by whisker program.
% Written by Erik Jacobsson 2022-03-09
% -----------------------------------

close all; clc; clear


%% Extracat and plot times

% Load times
time_dict = load_times();

% Print times
T = array2table(time_dict.values, 'RowNames', time_dict.keys, ...
    'VariableNames',{'Times'});
disp(T)
T2 = array2table([time_dict.keys, time_dict.values]);
writetable(T2, 'myData.csv', 'WriteRowNames', true);

%% Functions


function time_dict = load_times()
s = what('../single_study/mat_files');
file_location = s.path;
filename = 'times.mat';
filename = fullfile(file_location,filename);
load(filename,'exec_time', 'init_time','stiffness_time', ...
              'assem_stiffness_time','spaputval_time','solveq_time',...
              'force_time','assem_force_time','mat_time','plotvtk_time')

names = ["exec_time", "init_time", "stiffness_time", ...
         "assem_stiffness_time", "spaputval_time", "solveq_time", ...
         "force_time", "assem_force_time", "mat_time", "plotvtk_time"]';

times = [exec_time, init_time,stiffness_time, assem_stiffness_time,...
          spaputval_time, solveq_time, force_time, assem_force_time, ...
          mat_time, plotvtk_time]'/60;
time_dict = dictionary(names,times);

end