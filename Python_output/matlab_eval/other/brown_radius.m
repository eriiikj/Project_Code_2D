% Calfem
addpath /home/er7128ja/Nextcloud/Projekt/Matlab_Calfem/Calfem/fem
close all; clc; clear

% Geometry
width       = 2.5;
ngrains     = 2;
grain_width = width/ngrains;

% IMC steps
IMC_steps = 80;

% Define volume increase
V_IMC_over_A = 0.2;
vol_IMC      = 0;
radius_IMC   = 0;
vol_IMC_inc  = V_IMC_over_A/IMC_steps;

% Define plot variables
vol_inc_steps    = zeros(IMC_steps,1);
radius_inc_steps = zeros(IMC_steps,1);

% Increase volume
for k=1:IMC_steps

    % Update volume incrementally
    vol_IMC = vol_IMC + vol_IMC_inc;

    % Find corresponding IMC radius
    radius_IMC = sqrt(grain_width*2*vol_IMC/pi);

    % Save for plotting
    vol_inc_steps(k)    = vol_IMC;
    radius_inc_steps(k) = radius_IMC;

    disp('----')
    disp(['Load step ', num2str(k)])
    disp(['IMC volume per area: ', num2str(vol_IMC)])
    disp(['IMC radius: ', num2str(radius_IMC)])
end

% Plot radius versus IMC volumer over area
figure()
plot(radius_inc_steps, vol_inc_steps, 'ko-')
xlabel('r [um]')
ylabel('V_{imc}/A [um]')
title('IMC volume over area (r)')

% IMC growth over time
figure()
time_vec = (vol_inc_steps/0.014).^(1/0.6);
time_vec2 = (vol_inc_steps/0.014).^(1/1);
plot(time_vec, vol_inc_steps, 'ko-')
xlabel('t [hr]')
ylabel('V_{imc}/A [um]')
title('IMC volume over area (t)')

% % Plot radius versus IMC volumer over area
% plot(radius_vec, V_IMC_over_A_vec, 'r-', 'LineWidth',2)
% xlabel('r_{imc}')
% ylabel('V_{imc}/A [um]')
% title('IMC volume over area (r)')



% 
% % Radius vs time
% figure()
% plot(radius_vec, time_vec, 'b-', 'LineWidth',2)
% xlabel('r [um]')
% ylabel('t [hr]')
% title('IMC radius vs time')
