close all; clear; clc;


% Boundaries
ex = [0 2 2 0]; ey = [0 0 2 2];
patch(ex,ey,'white')
hold on

% Line coordinates
line_coord = [0.0 1.6; 0.1 1.43; 0.2 1.52; 0.9 1.44; 1.4 1.0; 2.0 0.9];

% Line conn OBS not sorted
line_conn = [1 2; 3 4; 2 3; 4 5; 5 6];

% Line ex and line ey
line_ex = reshape(line_coord(reshape(line_conn', 1, []),1),2,[])';
line_ey = reshape(line_coord(reshape(line_conn', 1, []),2),2,[])';

line_ex_o = line_ex;
line_ey_o = line_ey;

% Plot line segments
figure(1)
plot(line_ex',line_ey','-ko','LineWidth',2)


% Remove lines shorter than threhhold
lseg      = 0.3;  % Line segment threshhold
nseg      = size(line_ex,1);
rmseg     = zeros(nseg,1); % Info on if line seg is to be removed
rmseg     = rmseg>1;
rm        = true;
isegs     = (1:nseg)';
while (rm)
    rm = false;
    segloop = isegs(~rmseg);
    nloop   = sum(~rmseg);
    for k=1:nloop

        % Segment
        iseg = segloop(k);
    
        % 1) Take out vertex coordinates
        P1 = [line_ex_o(iseg,1),line_ey_o(iseg,1)];
        P2 = [line_ex_o(iseg,2),line_ey_o(iseg,2)];
    
        % 2) Compute line length
        line_length = norm(P1-P2);
    
        % 3) Check if line_length smaller than threshhold
        if (line_length<lseg)
    
            % Continue removing lines
            rm = true;

            % Check if P1 or P2 on boundary
            P1B = abs(P1(1)-min(ex))<1e-12 | abs(P1(1)-max(ex))<1e-12 | ...
                  abs(P1(2)-min(ey))<1e-12 | abs(P1(2)-max(ey))<1e-12;
            P2B = abs(P2(1)-min(ex))<1e-12 | abs(P2(1)-max(ex))<1e-12 | ...
                  abs(P2(2)-min(ey))<1e-12 | abs(P2(2)-max(ey))<1e-12;

            if (~P1B && ~P2B)
                % Neither P1 or P2 on domain boundary
                idx = abs(line_ex_o-P1(1))<1e-12 & abs(line_ey_o-P1(2))<1e-12 | ...
                      abs(line_ex_o-P2(1))<1e-12 & abs(line_ey_o-P2(2))<1e-12;            
                N   = (P1+P2)/2;
            elseif (P1B)
                % P1 on domain boundary. Set P2 at P1 pos.
                idx = abs(line_ex_o-P2(1))<1e-12 & abs(line_ey_o-P2(2))<1e-12;
                N   = P1;
            elseif(P2B)
                % P2 on domain boundary. Set P1 at P2 pos.
                idx = abs(line_ex_o-P1(1))<1e-12 & abs(line_ey_o-P1(2))<1e-12;
                N   = P2;
            end      
    
            % Update line coordinates
            line_ex(idx) = N(1);
            line_ey(idx) = N(2);          
        
        end    
    end
    % Extract lines
    rmseg     = rmseg | abs(line_ex(:,1)-line_ex(:,2))<1e-12 & abs(line_ey(:,1)-line_ey(:,2))<1e-12;
    nseg      = sum(~rmseg);
    line_ex_o = line_ex;
    line_ey_o = line_ey;
end

line_ex = line_ex(rmseg==0,:);
line_ey = line_ey(rmseg==0,:);



% Plot new line segments
figure(2)
patch(ex,ey,'white')
hold on
plot(line_ex',line_ey','-ko','LineWidth',2)
