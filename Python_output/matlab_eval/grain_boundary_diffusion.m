% --- Diffusion coefficient D for gb diffusion
close all; clc; clear;

% Bulk and gb diffusion coefficients
G_v  = 1e-9;
G_gb = G_v*10;


% Distance to gb
phi = linspace(-40e-5,40e-5);

% Coefficient
alpha = 1e4;

% G(phi)
G = G_v + (G_gb - G_v)*exp(-alpha*abs(phi));



% Plot
figure(1)
plot(phi, G, 'LineWidth',2)




