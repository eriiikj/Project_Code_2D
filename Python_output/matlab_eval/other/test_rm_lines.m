close all; clear; clc;


% Boundaries
ex = [0 2 2 0]; ey = [0 0 2 2];
patch(ex,ey,'white')
hold on

% Line coordinates
line_coord = [0.0 1.6; 0.3 1.1; 0.4 1.2; 0.47 1.2; 1.4 1.4; 2.0 0.9];

% Line conn OBS not sorted
line_conn = [1 3; 6 5; 2 3; 4 2; 4 5]

% Line ex and line ey
line_ex = reshape(line_coord(reshape(line_conn', 1, []),1),2,[])'
line_ey = reshape(line_coord(reshape(line_conn', 1, []),2),2,[])'

% Plot line segments
figure(1)
plot(line_ex',line_ey','-ko','LineWidth',2)



% Remove lines shorter than threhhold
lseg      = 0.2;  % Line segment threshhold
nseg      = size(line_ex,1);
nsegcoord = size(line_coord,1);
rm_seg    = zeros(nseg,1); % Info on if line seg is to be removed
rm_coord  = zeros(nsegcoord,1); % Info on if line coord is to be removed
for iseg=1:nseg

    % 1) Take out vertex indecies
    P1idx = line_conn(iseg,1);
    P2idx = line_conn(iseg,2);

    % 2) Take out vertex coordinates
    P1 = line_coord(P1idx,:);
    P2 = line_coord(P2idx,:);

    % 3) Compute line length
    line_length = norm(P1-P2);

    % 4) Check if line_length smaller than threshhold
    if (line_length<lseg)
        rm_seg(iseg)    = 1;  % Indicate that line is to be removed
        rm_coord(P2idx) = 1;  % Indicate that linecoord is to be removed

        % Update line_coord of P1 to N
        line_coord(P1idx,:) = (P1+P2)/2;

        % Change all P2idx to P1idx
        line_conn(line_conn==P2idx) = P1idx;

    end
end

% Remove line segments
% line_conn = line_conn(rm_seg==0,:);

% % Remove vertices
% line_coord    = line_coord(rm_coord==0,:);
% 

% 
% % Restructure line_conn to be compatible with new line coord
% ncoord          = sum(rm_coord==0);  % No. of vertices in new line_coord
% linecoord_idx   = find(rm_coord==0); % Old idx to be replaced
% line_coordidx_r = (1:ncoord)';       % New idx
% 
% % Iterate through each value in linecoord_idx and replace with corr. 
% % line_coordidx_r value
% for i = 1:length(linecoord_idx)
%     replace_idx            = (line_conn == linecoord_idx(i));
%     line_conn(replace_idx) = line_coordidx_r(i);
% end





% Line ex and line ey
line_ex = reshape(line_coord(reshape(line_conn', 1, []),1),2,[])';
line_ey = reshape(line_coord(reshape(line_conn', 1, []),2),2,[])';


% Plot line segments
figure(2)
patch(ex,ey,'white')
hold on
plot(line_ex',line_ey','-ko','LineWidth',2)