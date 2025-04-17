%% Problem 1

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2

%% Part a: Chief and deputy orbital elements

% For the chief
a_c = 36943; % km
e_c = 0.001;
% e_c = 0;
inc_c = deg2rad(59);
omega_c = deg2rad(188);
RAAN_c = deg2rad(84);
nu_c = deg2rad(0);
T = 2 * pi * sqrt(a_c^3 / mu);
n_c = sqrt(mu / a_c^3);

% For the deputy
delta_a = 0;
delta_e = - e_c * 0.0001 * a_c;
% delta_e = 0;
delta_i = deg2rad(0.1) * a_c;
% delta_i = 0;
delta_omega = 0;
delta_RAAN = - deg2rad(0.04) * a_c;
% delta_RAAN = 0;
delta_nu = deg2rad(-0.016) * a_c;
% delta_nu = 0;

%% Part b: Chief and deputy initial conditions

% Inertial position and velocity in ECI
a_d = a_c + delta_a;
e_d = e_c + delta_e / a_c;
inc_d = inc_c + delta_i / a_c;
omega_d = omega_c + delta_omega / a_c;
RAAN_d = RAAN_c + delta_RAAN / a_c;
nu_d = nu_c + delta_nu / a_c;

[pos_chief, vel_chief] = OEtoECI(a_c, e_c, inc_c, omega_c, RAAN_c, nu_c, mu);
[pos_deputy, vel_deputy] = OEtoECI(a_d, e_d, inc_d, omega_d, RAAN_d, nu_d, mu);

% Relative position and velocity in RTN
rho_ECI = pos_deputy - pos_chief;
rho_dot_ECI = vel_deputy - vel_chief;

rot_ECItoRTN = ECItoRTN([pos_chief; vel_chief]);
omega_vec = [0 0 sqrt(mu / (a_c^3 * (1 - e_c^2)^3)) * (1 + e_c * cos(nu_c))^2];

rho_RTN = rot_ECItoRTN * rho_ECI;
rho_dot_RTN = rot_ECItoRTN * rho_dot_ECI - cross(omega_vec', rho_RTN);

norm(rho_RTN) / (0.001 * a_c)

%% Part c: Integration constant

% Checking the condition for bounded motion
rho_dot_RTN(2) + 2 * n_c * rho_RTN(1)

% State Transition Matrix
mat_1 = [a_c * eye(3), zeros(3);
         zeros(3), a_c * n_c * eye(3)];
mat_2 = [1,    0,  1, 0, 0, 0;
         0,    2,  0, 1, 0, 0;
         0,    0,  0, 0, 0, 1;
         0,    1,  0, 0, 0, 0;
         -1.5, 0, -2, 0, 0, 0;
         0,    0,  0, 0, 1, 0];

% Computing the integration constant
% K = mat_2^(-1) * mat_1^(-1) * [rho_RTN / a_c; rho_dot_RTN / (a_c * n_c)];
K = mat_1^(-1) * [rho_RTN / a_c; rho_dot_RTN / (a_c * n_c)];

%% Part d: Propagation

N = 10000;
tspan = linspace(0, 15 * T, N);
relative_motion = zeros(N, 6);

for i=1:N
    relative_motion(i, :) = HCW_STM(n_c, tspan(i), K, mat_1);
end

figure
subplot(2,2,1)
plot(relative_motion(:, 2), relative_motion(:, 1))
grid on
axis equal
xlabel('T-axis')
ylabel('R-axis')

subplot(2,2,2)
plot(relative_motion(:, 3), relative_motion(:, 1))
grid on
axis equal
xlabel('N-axis')
ylabel('R-axis')

subplot(2,2,3)
plot(relative_motion(:, 2), relative_motion(:, 3))
grid on
axis equal
xlabel('T-axis')
ylabel('N-axis')

subplot(2,2,4)
plot3(relative_motion(:, 1), relative_motion(:, 2), relative_motion(:, 3))
grid on
axis equal
view(3)
xlabel('R-axis')
ylabel('T-axis')
zlabel('N-axis')


%% Functions

function state = HCW_STM(n, t, K, mat_1)
    s = sin(n * t);
    c = cos(n * t);
    mat_2 = [1, s, c, 0, 0, 0;
             -1.5 * n * t, 2 * c, -2 * s, 1, 0, 0;
             0, 0, 0, 0, s, c;
             0, c, -s, 0, 0, 0;
             -1.5, -2 * s, -2 * c, 0, 0, 0;
             0, 0, 0, 0, c, -s];
    state = mat_1 * mat_2 * K;
end