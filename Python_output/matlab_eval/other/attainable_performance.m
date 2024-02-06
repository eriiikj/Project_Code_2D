% -----------------------------------
% Attainable performance of c2dtl4_f
% Written by Erik Jacobsson 2022-12-06
% -----------------------------------
close all;clear;clc;

Bmem   = 205;    % Gbytes/s
Bfp_cc = 16;  % Flop/cc
Frek   = 2.25;   % Gcc/s
Bfp    = Bfp_cc*Frek;

Ifp  = 1472;
Imem = 3224;

AI   = Ifp/Imem;
bfpu = min(36,AI*Bmem);

AI_point = Bfp/Bmem;

AIvec1  = linspace(0,AI_point,100);
AIvec2  = linspace(AI_point,1,100);
bfpvec1 = AIvec1*Bmem;
bfpvec2 = Bfp*ones(100,1);

figure()
plot(AIvec1,bfpvec1,'k-',AIvec2,bfpvec2,'-k')
hold on
plot(AI,bfpu,'ko')
xlabel('AI [Flop/Byte]')
ylabel('Attainable performance [GFlop/s]')
axis([0 1 0 40])
box off


%% Time
nelm            = 14875;
nstep           = 10;
deltat_e        = Ifp/bfpu*1e-9;
delta_t         = nstep*nelm*deltat_e;
real_time       = 0.25;
attainable_perc = delta_t/real_time*100