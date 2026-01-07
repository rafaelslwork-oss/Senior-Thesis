clear; clc;

% Constants
G = 6.67430e-11;
Msun = 1.98847e30;
Mbh = 4.152e6 * Msun;
AU = 1.495978707e11;
pc = 3.085677581e16;

% Orbit parameters
e_DM = 0.88466;
p_AU = 420.73/2;
p = p_AU * AU;
L = sqrt(G*Mbh * p);

% Dark matter spike
Rs = 0.001 * pc;
rho0 = 1e5 * Msun / pc^3 ;
a = 2*pi*G*rho0*Rs;

% Initial conditions
r0 = p/(1 + e_DM);
vr0 = 0;
theta0 = 0;

% Equations of motion
odefun = @(t,y)[
    y(2); 
    (L^2./y(1).^3) - (G*Mbh./y(1).^2) - a; % dv_r/dt
    L./y(1).^2 
];

% Orbit integration
T_end = 2.0e9;
y0 = [r0; vr0; theta0];

Npts = 40000;
tspan = linspace(0, T_end, Npts);

options = odeset('RelTol',1e-12,'AbsTol',1e-12,'MaxStep',1e8);
[t_sol, y_sol] = ode45(odefun, tspan, y0, options);

r_DM = y_sol(:,1);
th_DM = y_sol(:,3);

r_DM_AU = r_DM/AU;
x_DM = r_DM_AU .* cos(th_DM);
y_DM = r_DM_AU .* sin(th_DM);

% Pure ellipse
e_ell = 0.8854;
theta_ell = linspace(0, 2*pi, 2000);
r_ell_AU = p_AU./(1 + e_ell.*cos(theta_ell));
x_ell = r_ell_AU .* cos(theta_ell);
y_ell = r_ell_AU .* sin(theta_ell);

% Pulling pure schwarzshcild from sageMath code
data = readtable("pureSchwarzschild.csv");
r_sage = data.r(:)/AU;
theta_sage = data.phi(:);
x_Sage = r_sage .* cos(theta_sage(:));
y_Sage = r_sage .* sin(theta_sage(:));

% Pulling Schwarzschild + NFW from sageMath code
data2 = readtable("SchwarzschildPlusNFW.csv");
r_sage2 = data2.r(:)/AU;
theta_sage2 = data2.phi(:);
x_Sage2 = r_sage2 .* cos(theta_sage2(:));
y_Sage2 = r_sage2 .* sin(theta_sage2(:));

% Plot
figure; hold on;
plot(x_ell,y_ell,'k','LineWidth',1);
plot(x_Sage2,y_Sage2,'m--','LineWidth',1);
plot(0,0,'gs','MarkerFaceColor','g','MarkerSize',6);

legend('Observed Orbit','Schwarzschild + NFW','Sgr A*','Location','southoutside');
xlabel('x (AU)');
ylabel('y (AU)');
title('S2 Orbit: Observed vs Schwarzschild + NFW');

axis equal; grid on;
xlim([-2500 500]);
ylim([-1000 1000]);

print(gcf,'InitialS2_Orbit.pdf','-dpdf','-bestfit');