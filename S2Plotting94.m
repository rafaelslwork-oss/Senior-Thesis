clear; clc;

% Physical constants
G   = 6.67430e-11;
c   = 299792458;
Msun = 1.98847e30;
Mbh  = 4.152e6 * Msun;
AU   = 1.495978707e11;
pc   = 3.085677581e16;

% Orbital elements for S2
e     = 0.88466;
p_AU  = 226.16;            % semilatus rectum giving rp = 120 AU
p     = p_AU * AU;
L     = sqrt(G*Mbh*p);

% Dark matter parameters
rho0 = 1e5 * Msun / pc^3; 
Rs   = 0.001 * pc;
a    = 2 * pi * G * rho0 * Rs;

% Initial conditions
r0   = p / (1 + e);
dr0  = 0;
theta0 = 0;

%Energy equations and calculation of equations of motion, relevant math
%done pen and paper
E_NoDM = 0.5*dr0^2 + L^2/(2*r0^2) - G*Mbh/r0;
E_DM   = 0.5*dr0^2 + L^2/(2*r0^2) - G*Mbh/r0 + a*r0;

P_NoDM = @(r) 2*E_NoDM.*r.^2 + 2*G*Mbh.*r - L^2;
P_DM   = @(r) -2*a.*r.^3 + 2*E_DM.*r.^2 + 2*G*Mbh.*r - L^2;

dt = 1e11;
N  = 8000;

r_DM  = zeros(1,N);  th_DM  = zeros(1,N);
r_NoDM= zeros(1,N);  th_NoDM= zeros(1,N);

r_DM(1)   = r0;  r_NoDM(1)   = r0;
th_DM(1)  = theta0;  th_NoDM(1) = theta0;

for k = 1:N-1
    rr1 = r_DM(k);
    rr2 = r_NoDM(k);

    % Radial derivatives
    dr_DM   = sqrt(P_DM(rr1)) / rr1;
    dr_NoDM = sqrt(P_NoDM(rr2)) / rr2;

    r_DM(k+1)   = rr1 + dr_DM*dt;
    r_NoDM(k+1) = rr2 + dr_NoDM*dt;

    th_DM(k+1)  = th_DM(k)  + (L/rr1^2)*dt;
    th_NoDM(k+1)= th_NoDM(k)+ (L/rr2^2)*dt;
end

% Convert to AU
r_DM_AU   = r_DM / AU;
r_NoDM_AU = r_NoDM / AU;

x_DM   = r_DM_AU   .* cos(th_DM);
y_DM   = r_DM_AU   .* sin(th_DM);
x_NoDM = r_NoDM_AU .* cos(th_NoDM);
y_NoDM = r_NoDM_AU .* sin(th_NoDM);

%Schwarzschild from sageMath
data = readtable("pureSchwarzschild.csv");
r_sage    = data.r(:) / AU;
theta_sage = data.phi(:);
x_Sage = r_sage .* cos(theta_sage(:));
y_Sage = r_sage .* sin(theta_sage(:));

%Plotting
figure; hold on;
plot(x_NoDM,  y_NoDM,  'k--', 'LineWidth', 2);     
plot(x_DM,    y_DM,    'r',   'LineWidth', 1.5);   
plot(x_Sage,  y_Sage,  'b',   'LineWidth', 1.5);   
plot(0,0,'gs','MarkerFaceColor','g','MarkerSize',6);

xlabel('x (AU)');
ylabel('y (AU)');
title('Three-Model Orbit Comparison for S2: Newtonian / Newtonian+DM / GR');
legend('Newtonian (No DM)', ...
       'Newtonian + DM Spike', ...
       'Schwarzschild GR (Sage)', ...
       'Sgr A*', 'Location','bestoutside');

axis equal; grid on;
xlim([-3000 2000]); ylim([-1000 1000]);

print(gcf,'S2_ThreeModelComparison.pdf','-dpdf','-bestfit');
