% --- Illustration of the interpolation function ---

close all; clear; clc;

% Signed distance function
a_ls1 = linspace(-10,10)';
a_ls2 = linspace(10,-10)';

% Interpolation function
w    = 5;
phi1 = 0.5*(1 - tanh((a_ls1)/(w/2)));
phi2 = 0.5*(1 - tanh((a_ls2)/(w/2)));
h1   = phi1.^2./(phi1.^2 + phi2.^2);
h2   = phi2.^2./(phi1.^2 + phi2.^2);
h11  = a_ls1.^2./(a_ls1.^2 + a_ls2.^2);
h22  = a_ls2.^2./(a_ls1.^2 + a_ls2.^2);
% phi1b = 0.5*(1 - tanh((a_ls1+w*0.5)/(w/4/2)));

% Plot interpolation function
figure(1)
plot(a_ls1, a_ls1, 'LineWidth',2,'DisplayName','als1')
hold on
plot(a_ls1, a_ls2, 'LineWidth',2,'DisplayName','als2')
plot(a_ls1, phi1, 'LineWidth',2,'DisplayName','phi1')

plot(a_ls1, phi2, 'LineWidth',2,'DisplayName','phi2')
plot(a_ls1, h1, 'LineWidth',2,'DisplayName','h1')
plot(a_ls1, h11, 'LineWidth',2,'DisplayName','h11')
ylabel('Interpolation function')
xlabel('Signed distance')
title('Interpolation function')
legend
axis([-10 10 0 1])

% % Interpolation function2
% w    = 5;
% phi1 = 0.5*(1 - tanh((a_ls1+w*0.5)/(w/2)));
% phi2 = 0.5*(1 - tanh((a_ls2+w*0.5)/(w/2)));
% h1   = phi1.^2./(phi1.^2 + phi2.^2);
% 
% % Plot interpolation function
% figure(2)
% plot(a_ls1, phi1, 'LineWidth',2,'DisplayName','phi1')
% hold on
% plot(a_ls1, phi2, 'LineWidth',2,'DisplayName','phi2')
% plot(a_ls1, h1, 'LineWidth',2,'DisplayName','h1')
% ylabel('Interpolation function')
% xlabel('Signed distance')
% title('Interpolation function')
% legend


% % --- Example ---
% 
% % Stresses
% sigma1 = 50;
% sigma2 = 100; 
% sigma3 = 50;
% 
% sigma_max = max([sigma1,sigma2,sigma3])
% 
% % Signed distances to point
% a_ls1 = 1;
% a_ls2 = -1;
% a_ls3 = 2;
% 
% % Corresponding phi's
% phi1 =  0.5*(1 - tanh(a_ls1/(w/2)))
% phi2 =  0.5*(1 - tanh(a_ls2/(w/2)))
% phi3 =  0.5*(1 - tanh(a_ls3/(w/2)))
% 
% % Corresponding h's (compare JH)
% h1 = phi1^2/(phi1^2 + phi2^2 + phi3^2)
% h2 = phi2^2/(phi1^2 + phi2^2 + phi3^2)
% h3 = phi3^2/(phi1^2 + phi2^2 + phi3^2)
% 
% h1 = phi1/(phi1 + phi2 + phi3)
% h2 = phi2/(phi1 + phi2 + phi3)
% h3 = phi3/(phi1 + phi2 + phi3)
% 
% 
% phi_sum = phi1 + phi2 + phi3
% h_sum   = h1 + h2 + h3
% 
% sigma = h1*sigma1 + h2*sigma2 + h3*sigma3
