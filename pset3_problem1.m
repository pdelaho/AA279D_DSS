%% Problem 1

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2

%% Part a: Chief and deputy orbital elements

% For the chief
a_c = 36943; % km
e_c = 0.0001;
inc_c = deg2rad(59);
omega_c = deg2rad(188);
RAAN_c = deg2rad(84);
nu_c = deg2rad(0);
T = 2 * pi * sqrt(a_c^3 / mu);
n_c = sqrt(mu / a_c^3);

% For the deputy
delta_a = 0;
delta_e = - e_c * 0.1 * a_c;
delta_i = deg2rad(0.0005) * a_c;
delta_omega = 0;
delta_RAAN = - deg2rad(0.0001) * a_c;
delta_nu = deg2rad(-0.0016) * a_c;

%% Part b: Chief and deputy initial conditions

% Inertial position and velocity in ECI
a_d = a_c + delta_a;
e_d = e_c + delta_e / a_c;
inc_d = inc_c + delta_i / a_c;
omega_d = omega_c + delta_omega / a_c;
RAAN_d = RAAN_c + delta_RAAN / a_c;
nu_d = nu_c + delta_nu / a_c;

[pos_chief, vel_chief] = OE2ECI(a_c, e_c, inc_c, omega_c, RAAN_c, nu_c, mu);
[pos_deputy, vel_deputy] = OE2ECI(a_d, e_d, inc_d, omega_d, RAAN_d, nu_d, mu);

% Relative position and velocity in ECI
rho_ECI = pos_deputy - pos_chief;
rho_dot_ECI = vel_deputy - vel_chief;

% Relative position and velocity in RTN
state_RTN = ECI2RTN([pos_chief; vel_chief], [rho_ECI; rho_dot_ECI], mu);
rho_RTN = state_RTN(1:3)';
rho_dot_RTN = state_RTN(4:6)';

norm(rho_RTN) / (a_c)

%% Part c: Integration constant

% State Transition Matrix at t=0
scaling = [a_c * eye(3), zeros(3);
         zeros(3), a_c * n_c * eye(3)];

STM = [   1, 0,  1, 0, 0, 0;
          0, 2,  0, 1, 0, 0;
          0, 0,  0, 0, 0, 1;
          0, 1,  0, 0, 0, 0;
       -1.5, 0, -2, 0, 0, 0;
          0, 0,  0, 0, 1, 0];

% Computing the integration constant
K = STM^(-1) * scaling^(-1) * [rho_RTN; rho_dot_RTN]

%% Part d: Propagation

N = 15000;
tspan = linspace(0, 15 * T, N);
relative_motion = zeros(N, 6);

for i=1:N
    STM_temporary = HCW_STM(n_c, tspan(i));
    relative_motion(i, :) = scaling * STM_temporary * K;
end

% Shape of the ellipses
a_c * K(1)
a_c * K(4)
a_c * sqrt(K(2)^2 + K(3)^2)
2 * a_c * sqrt(K(2)^2 + K(3)^2)
a_c * sqrt(K(5)^2 + K(6)^2)

figure
subplot(2,2,1)
plot(relative_motion(:, 2), relative_motion(:, 1))
grid on
axis equal
xlabel('T-axis [km]')
ylabel('R-axis [km]')

subplot(2,2,2)
plot(relative_motion(:, 3), relative_motion(:, 1))
grid on
axis equal
xlabel('N-axis [km]')
ylabel('R-axis [km]')

subplot(2,2,3)
plot(relative_motion(:, 2), relative_motion(:, 3))
grid on
axis equal
xlabel('T-axis [km]')
ylabel('N-axis [km]')

subplot(2,2,4)
plot3(relative_motion(:, 1), relative_motion(:, 2), relative_motion(:, 3))
grid on
axis equal
view(3)
xlabel('R-axis [km]')
ylabel('T-axis [km]')
zlabel('N-axis [km]')

figure
subplot(2,2,1)
plot(relative_motion(:, 5), relative_motion(:, 4))
grid on
axis equal
xlabel('T-axis [km/s]')
ylabel('R-axis [km/s]')

subplot(2,2,2)
plot(relative_motion(:, 6), relative_motion(:, 4))
grid on
axis equal
xlabel('N-axis [km/s]')
ylabel('R-axis [km/s]')

subplot(2,2,3)
plot(relative_motion(:, 5), relative_motion(:, 6))
grid on
axis equal
xlabel('T-axis [km/s]')
ylabel('N-axis [km/s]')

subplot(2,2,4)
plot3(relative_motion(:, 4), relative_motion(:, 5), relative_motion(:, 6))
grid on
axis equal
view(3)
xlabel('R-axis [km/s]')
ylabel('T-axis [km/s]')
zlabel('N-axis [km/s]')

figure
plot(relative_motion(:, 2), relative_motion(:, 1))
grid on
axis equal
xlim([-1.82, -1.74])
xlabel('T-axis [km]')
ylabel('R-axis [km]')

figure
plot(relative_motion(:, 2), relative_motion(:, 1))
grid on
axis equal
xlabel('T-axis [km]')
ylabel('R-axis [km]')

figure
plot(relative_motion(:, 3), relative_motion(:, 1))
grid on
axis equal
xlabel('N-axis [km]')
ylabel('R-axis [km]')

%% Part e: 

% Checking the condition for stable motion
rho_dot_RTN(2) + 2 * n_c * rho_RTN(1)

% Enforcing the stable motion in the along-track direction
rho_RTN(1) = - rho_dot_RTN(2) / (2 * n_c);

% Computing the new integration constants
K_bis = STM^(-1) * scaling^(-1) * [rho_RTN; rho_dot_RTN]

N = 15000;
tspan = linspace(0, 15 * T, N);
relative_motion_bis = zeros(N, 6);

for i=1:N
    STM_temporary = HCW_STM(n_c, tspan(i));
    relative_motion_bis(i, :) = scaling * STM_temporary * K_bis;
end

figure
subplot(2,2,1)
plot(relative_motion_bis(:, 2), relative_motion_bis(:, 1))
grid on
axis equal
xlabel('T-axis [km]')
ylabel('R-axis [km]')

subplot(2,2,2)
plot(relative_motion_bis(:, 3), relative_motion_bis(:, 1))
grid on
axis equal
xlabel('N-axis [km]')
ylabel('R-axis [km]')

subplot(2,2,3)
plot(relative_motion_bis(:, 2), relative_motion_bis(:, 3))
grid on
axis equal
xlabel('T-axis [km]')
ylabel('N-axis [km]')

subplot(2,2,4)
plot3(relative_motion_bis(:, 1), relative_motion_bis(:, 2), relative_motion_bis(:, 3))
grid on
axis equal
view(3)
xlabel('R-axis [km]')
ylabel('T-axis [km]')
zlabel('N-axis [km]')

figure
subplot(2,2,1)
plot(relative_motion_bis(:, 5), relative_motion_bis(:, 4))
grid on
axis equal
xlabel('T-axis [km/s]')
ylabel('R-axis [km/s]')

subplot(2,2,2)
plot(relative_motion_bis(:, 6), relative_motion_bis(:, 4))
grid on
axis equal
xlabel('N-axis [km/s]')
ylabel('R-axis [km/s]')

subplot(2,2,3)
plot(relative_motion_bis(:, 5), relative_motion_bis(:, 6))
grid on
axis equal
xlabel('T-axis [km/s]')
ylabel('N-axis [km/s]')

subplot(2,2,4)
plot3(relative_motion_bis(:, 4), relative_motion_bis(:, 5), relative_motion_bis(:, 6))
grid on
axis equal
view(3)
xlabel('R-axis [km/s]')
ylabel('T-axis [km/s]')
zlabel('N-axis [km/s]')

%% Functions

function M = HCW_STM(n, t)
    s = sin(n * t);
    c = cos(n * t);
    M = [           1,      s,      c, 0, 0,  0;
         -3/2 * n * t,  2 * c, -2 * s, 1, 0,  0;
                    0,      0,      0, 0, s,  c;
                    0,      c,     -s, 0, 0,  0;
                 -3/2, -2 * s, -2 * c, 0, 0,  0;
                    0,      0,      0, 0, c, -s];
end