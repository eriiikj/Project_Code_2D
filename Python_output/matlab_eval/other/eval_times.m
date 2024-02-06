% -----------------------------------
% Eval times from json
% Written by Erik Jacobsson 2022-12-02
% -----------------------------------

close all; clear; clc;


% Change to file location
file_path = ['/home/er7128ja/Nextcloud/Kurser/PDC_summer_school'...
             '/Project/Project_Code/Python_output/single_study'];

FileList = dir(fullfile(file_path,'times*'));
nfiles = length(FileList);

% Initiate nfiles
fname = fullfile(file_path, FileList(1).name); 
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);
names = fieldnames(val);
ntimes = length(fieldnames(val));

xlim = 22;
ylim = 14;

%% Execution times

exec_times_threads = zeros(nfiles,2);

for j=1:nfiles
    fname = fullfile(file_path, FileList(j).name);
    cores = str2double(FileList(j).name(7:end-5));
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);
    exec_times_threads(j,1) = cores;
    exec_times_threads(j,2) = val.exec_time;
end
[~,idx] = sort(exec_times_threads(:,1));
exec_times_threads = exec_times_threads(idx,:);
exec_speedup       = exec_times_threads(1,2)./ ...
                     exec_times_threads(:,2);

% % Plot
% figure()
% plot(exec_times_threads(:,1), exec_speedup, ...
%      'ko', 'MarkerFaceColor',[0.5, 0.5, 0.5])
% axis([0 xlim 0 ylim])
% box off
% legend('Measured speedup', 'Location','best')
% xlabel('Number of threads')
% xticks(2:2:xlim-1)
% ylabel('Speedup')
% title('Execution Time Speedup')

%% Stiffness times

stiffness_times_threads = zeros(nfiles,2);

for j=1:nfiles
    fname = fullfile(file_path, FileList(j).name);
    cores = str2double(FileList(j).name(7:end-5));
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);
    stiffness_times_threads(j,1) = cores;
    stiffness_times_threads(j,2) = val.stiffness_time;
end
[~,idx] = sort(stiffness_times_threads(:,1));
stiffness_times_threads = stiffness_times_threads(idx,:);
stiffness_speedup       = stiffness_times_threads(1,2)./ ...
                          stiffness_times_threads(:,2);

% Plot
figure()
plot(stiffness_times_threads(:,1), stiffness_speedup, ...
     'ko', 'MarkerFaceColor',[0.5, 0.5, 0.5])
axis([0 xlim 0 ylim])
box off
legend('Measured speedup', 'Location','best')
xlabel('Number of threads')
xticks(2:2:xlim-1)
ylabel('Speedup')
title('Stiffness Speedup')

%% Assem stiffness times

assem_stiffness_times_threads = zeros(nfiles,2);

for j=1:nfiles
    fname = fullfile(file_path, FileList(j).name);
    cores = str2double(FileList(j).name(7:end-5));
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);
    assem_stiffness_times_threads(j,1) = cores;
    assem_stiffness_times_threads(j,2) = val.assem_stiffness_time;
end
[~,idx] = sort(assem_stiffness_times_threads(:,1));
assem_stiffness_times_threads = assem_stiffness_times_threads(idx,:);
assem_stiffness_speedup       = assem_stiffness_times_threads(1,2)./ ...
                                assem_stiffness_times_threads(:,2);

% Plot
figure()
plot(assem_stiffness_times_threads(:,1), assem_stiffness_speedup,...
     'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5])
axis([0 xlim 0 ylim])
box off
legend('Measured speedup', 'Location','best')
xlabel('Number of threads')
xticks(2:2:xlim-1)
ylabel('Speedup')
title('Assem Stiffness Speedup')

%% Solveq times

solveq_times_threads = zeros(nfiles,2);

for j=1:nfiles
    fname = fullfile(file_path, FileList(j).name);
    cores = str2double(FileList(j).name(7:end-5));
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);
    solveq_times_threads(j,1) = cores;
    solveq_times_threads(j,2) = val.solveq_time;
end
[~,idx] = sort(solveq_times_threads(:,1));
solveq_times_threads = solveq_times_threads(idx,:);
solveq_speedup       = solveq_times_threads(1,2)./ ...
                       solveq_times_threads(:,2);

% Plot
figure()
plot(solveq_times_threads(:,1), solveq_speedup,...
     'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5])
axis([0 xlim 0 ylim])
box off
legend('Measured speedup', 'Location','best')
xlabel('Number of threads')
xticks(2:2:xlim-1)
ylabel('Speedup')
title('Solveq Speedup')



%% Force times

force_times_threads = zeros(nfiles,2);

for j=1:nfiles
    fname = fullfile(file_path, FileList(j).name);
    cores = str2double(FileList(j).name(7:end-5));
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);
    force_times_threads(j,1) = cores;
    force_times_threads(j,2) = val.force_time;
end
[~,idx] = sort(force_times_threads(:,1));
force_times_threads = force_times_threads(idx,:);
force_speedup       = force_times_threads(1,2)./ ...
                      force_times_threads(:,2);

% Plot
figure()
plot(force_times_threads(:,1), force_speedup, ...
     'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5])
axis([0 xlim 0 ylim])
box off
legend('Measured speedup', 'Location','best')
xlabel('Number of threads')
xticks(2:2:xlim-1)
ylabel('Speedup')
title('Force Speedup')


%% Assem force times

assem_force_times_threads = zeros(nfiles,1);

for j=1:nfiles
    fname = fullfile(file_path, FileList(j).name);
    cores = str2double(FileList(j).name(7:end-5));
    fid = fopen(fname); 
    raw = fread(fid,inf); 
    str = char(raw'); 
    fclose(fid); 
    val = jsondecode(str);
    assem_force_times_threads(j,1) = cores;
    assem_force_times_threads(j,2) = val.assem_force_time;
end
[~,idx] = sort(assem_force_times_threads(:,1));
assem_force_times_threads = assem_force_times_threads(idx,:);
assem_force_speedup       = assem_force_times_threads(1,2)./ ...
                            assem_force_times_threads(:,2);

% Plot
figure()
plot(assem_force_times_threads(:,1), assem_force_speedup, ...
     'ko', 'MarkerFaceColor', [0.5, 0.5, 0.5])
axis([0 xlim 0 ylim])
box off
legend('Measured speedup', 'Location','best')
xlabel('Number of threads')
xticks(2:2:xlim-1)
ylabel('Speedup')
title('Assem Force Speedup')