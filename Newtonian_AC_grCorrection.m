% This code is plotting the orbit of S2 with Newtonian Mechanics + an
% Adiabatically Compressed Dark Matter Halo Model (Gondolo & Silk) + a
% first order Post-Newtonian general relativity correction term. The
% analytical math was done pen and paper.

clear; clc;
% Constants
G = 6.67430e-11;
Msun = 1.98847e30;
Mbh = 4.30e6 * Msun;
AU = 1.495978707e11;
pc = 3.085677581e16;
c = 299792458;

% Values and plotting for observed orbit
e_DM = 0.8854;
r_peri_AU = 120;

a_AU = r_peri_AU/(1 - e_DM);
a = a_AU * AU;
p_AU = r_peri_AU*(1 + e_DM);
p = p_AU * AU;

% DM Spike parameters
gamma_sp = 7/3;
rho_sp = 1e-2 * Msun / pc^3 * 7.4e12;
r_sp = .5 * pc;
A = (4*pi*rho_sp*r_sp^gamma_sp)/(3 - gamma_sp);

% Angular momentum calculation
r0 = r_peri_AU * AU;
vr0 = 0;
theta0 = 0;

v0 = (1 + e_DM) * sqrt(G*Mbh/p);
L  = r0 * v0;

% Equations of motion calculation
odefun = @(t,y)[ ...
    y(2); ...
    (L^2./y(1).^3) - (G*Mbh./y(1).^2) - (G*A.*y(1).^(1 - gamma_sp)) - (3*G*Mbh*L^2./(c^2.*y(1).^4)); ...
    L./y(1).^2 ...
];

T_end = 360 * 15.9 * 365.25*24*3600; % integration end time
y0 = [r0; vr0; theta0];

Npts = 1e5;
tspan = linspace(0, T_end, Npts);

options = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',1e8);
[t_sol, y_sol] = ode45(odefun, tspan, y0, options);

% Conversion to Cartesians
r_DM = y_sol(:,1);
th_DM = y_sol(:,3);

r_DM_AU = r_DM/AU;
x_DM = r_DM_AU.*cos(th_DM);
y_DM = r_DM_AU.*sin(th_DM);

% e_ell = 0.8831;
% theta_ell = linspace(0, 2*pi, 5e6);
% r_ell_AU = p_AU./(1 + e_ell.*cos(theta_ell));
% x_ell = r_ell_AU.*cos(theta_ell);
% y_ell = r_ell_AU.*sin(theta_ell);

figure; hold on;
%plot(x_ell,y_ell,'k','LineWidth',1.5);
plot(x_DM,y_DM,'r','LineWidth', .05);
plot(0,0,'gs','MarkerFaceColor','g','MarkerSize',6);

%Plotting Observed Orbit
%data4 = readmatrix('ObservedOrbit.csv');
% x4 = data4(:,1);
% y4 = data4(:,2);

%plot(x4, y4, 'g', 'LineWidth', 1);

legend('Newtonian + GR Correction + AC DM Halo', 'Sgr A*', 'Location','southoutside');

xlabel('x (AU)');
ylabel('y (AU)');
title('S2 Orbit: Newtonian + GR Correction Precession', 'FontSize', 22);

axis equal; grid on;
xlim([-2000 300]);
ylim([-1000 600]);

hold off

% Saving for future plotting
print(gcf,'Newtonian_GR_DM_Orbit','-dpdf','-bestfit');
output = [x_DM(:), y_DM(:)];
writematrix(output,'Newtonian_GR_DM_Orbit.csv');

% Precession Calculation
r   = y_sol(:,1);
phi = unwrap(y_sol(:,3));

peri_idx = find(r(2:end-1) < r(1:end-2) & r(2:end-1) < r(3:end)) + 1;

phi_peri = phi(peri_idx);
dphi = diff(phi_peri); 

dphi_last = dphi(end);
delta_rad = dphi_last - 2*pi; 

delta_deg = delta_rad * 180/pi;

% RMS sigma calculation
theta_DM_unwrapped = unwrap(th_DM);
theta_DM_unwrapped = theta_DM_unwrapped(:);
r_DM_AU_vec = r_DM_AU(:);

r_DM_interp = interp1(theta_DM_unwrapped, r_DM_AU_vec, theta_ell(:), 'linear', 'extrap');

x_DM_interp = r_DM_interp .* cos(theta_ell(:));
y_DM_interp = r_DM_interp .* sin(theta_ell(:));

x_ell_vec = x_ell(:);
y_ell_vec = y_ell(:);

dist = sqrt((x_DM_interp - x_ell_vec).^2 + (y_DM_interp - y_ell_vec).^2);

max_diff = max(dist);
sigma = sqrt(mean(dist.^2));
