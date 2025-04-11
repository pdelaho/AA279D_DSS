%% Problem 1  

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2

%% Part a: Chief and deputy orbits

% For the chief
a_c = 36943; % km, semi-major axis
e_c = 0.8111; % eccentricity
inc_c = deg2rad(59); % inclination
omega_c = deg2rad(188); % argument of periapsis
RAAN_c = deg2rad(84); % right ascension of the ascending node
nu_c = 0; % true anomaly
T = 2 * pi * sqrt(a_c^3 / mu);

% For the deputy
% Write it in terms of delta in [km]
delta_a = 0;
delta_e = 0;
delta_i = deg2rad(0.1) * a_c;
delta_omega = deg2rad(0.05) * a_c;
delta_RAAN = 0;
delta_nu = deg2rad(-0.016) * a_c;

a_d = a_c + delta_a;
e_d = e_c + delta_e / a_c;
inc_d = inc_c + delta_i / a_c;
omega_d = omega_c + delta_omega / a_c;
RAAN_d = RAAN_c + delta_RAAN / a_c;
nu_d = nu_c + delta_nu / a_c;

%% Part b: Numerical integration of NERM

[pos_c, vel_c] = OEtoECI(a_c, e_c, inc_c, omega_c, RAAN_c, nu_c, mu);
[pos_d, vel_d] = OEtoECI(a_c, e_d, inc_d, omega_d, RAAN_d, nu_d, mu);

init_cond = [pos_c, vel_c, pos_d - pos_c, vel_d - vel_c];
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 100);
N = 10000;
tspan = linspace(0, 10 * T, N);
[t, y] = ode89(@(t, state) NERM(t, state, mu), tspan, init_cond, options);

relative_RTN = zeros(N, 6);

for i=1:N
    pos_c = y(i, 1:3);
    vel_c = y(i, 4:6);
    rho = y(i, 7:9);
    rho_dot = y(i, 10:12);
    [a, e, inc, omega, RAAN, nu] = Keplerian_elements(y(i, 1:6), mu);

    % Computing the RTN frame centered on the chief spacecraft
    R_vec = pos_c / norm(pos_c);
    N_vec = cross(pos_c, vel_c) / norm(cross(pos_c, vel_c));
    T_vec = cross(N_vec, R_vec);
    rot = [R_vec; T_vec; N_vec];
    theta_dot = [0, 0, - sqrt(mu / (a^3 * (1 - e^2)^3)) * (1 + e * cos(nu))^2];
    
    % Transforming the position and velocities in the RTN frame
    rho_RTN = rot * rho';
    rho_dot_RTN = rot * rho_dot' + cross(theta_dot, rho_RTN)';

    relative_RTN(i, 1:3) = rho_RTN;
    relative_RTN(i, 4:6) = rho_dot_RTN;
end

figure
subplot(3,2,1)
plot(t / 3600, relative_RTN(:, 1))
grid on
ylabel('R-component of position [km]')

subplot(3,2,3)
plot(t / 3600, relative_RTN(:, 2))
grid on
ylabel('T-component of position [km]')

subplot(3,2,5)
plot(t / 3600, relative_RTN(:, 3))
grid on
ylabel('N-component of position [km]')
xlabel('Time [hours]')

subplot(3,2,2)
plot(t / 3600, relative_RTN(:, 4))
grid on
ylabel('R-component of velocity [km/s]')

subplot(3,2,4)
plot(t / 3600, relative_RTN(:, 5))
grid on
ylabel('T-component of velocity [km/s]')

subplot(3,2,6)
plot(t / 3600, relative_RTN(:, 6))
grid on
ylabel('N-component of velocity [km/s]')
xlabel('Time [hours]')

figure
plot3(y(:, 7), y(:, 8), y(:, 9))
axis equal
grid on
xlabel('R-axis [km]')
ylabel('T-axis [km]')
zlabel('N-axis [km]')
title('Deputy orbit in RTN frame centered on chief')

%% Functions
function [pos_inertial, vel_inertial] = OEtoECI(a, e, inc, omega, RAAN, true_anom, mu)
    cosE = (e + cos(true_anom)) / (1 + e * cos(true_anom));
    sinE = sin(true_anom) * sqrt(1 - e^2) / (1 + e * cos(true_anom));
    n = sqrt(mu / a^3);
    
    % Computing the position and velocity vectors in the perifocal frame
    pos_perifocal = [a * (cosE - e) a * sqrt(1 - e^2) * sinE 0];
    vel_perifocal = a * n / (1 - e * cosE) * [-sinE sqrt(1 - e^2)*cosE 0];
    
    % Computing the rotation matrix between perifocal and inertial frame
    rotRAAN = [cos(RAAN) sin(-RAAN) 0;
               -sin(-RAAN) cos(RAAN) 0;
               0 0 1];
    roti = [1 0 0;
            0 cos(inc) sin(-inc);
            0 -sin(-inc) cos(inc)];
    rotomega = [cos(omega) sin(-omega) 0;
                -sin(-omega) cos(omega) 0;
                0 0 1];
    rot_perifocalTOinertial = rotRAAN * roti * rotomega;

    % Rotating the position and velocity vectors
    pos_inertial = rot_perifocalTOinertial * pos_perifocal';
    vel_inertial = rot_perifocalTOinertial * vel_perifocal';
end

function statedot = NERM(t, state, mu)
    statedot = zeros(size(state));
    statedot(1:3) = state(4:6);
    statedot(7:9) = state(10:12);

    r0 = state(1:3);
    rho = state(7:9);

    statedot(4:6) = - mu * r0 / norm(r0)^3;
    statedot(10:12) = - mu * (r0 + rho) / norm(r0 + rho)^3 + mu * r0 / norm(r0)^3;
end

function [a, e, i, omega, RAAN, nu] = Keplerian_elements(state, mu)
    pos = state(1:3);
    vel = state(4:6);
    h = cross(pos, vel);
    W = h / norm(h);
    i = atan2(sqrt(W(1)^2 + W(2)^2), W(3));
    RAAN = atan2(W(1), - W(2));
    p = norm(h)^2 / mu;
    a = 1 / (2 / norm(pos) - norm(vel)^2 / mu);
    n = sqrt(mu / a^3);
    e = sqrt(1 - p / a);
    E = atan2((dot(pos, vel) / (a^2 * n)), (1 - norm(pos)/a));
    nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
%     M = wrapTo2Pi(E - e * sin(E));
    u = atan2((pos(3) / sin(i)), (pos(1) * cos(RAAN) + pos(2) * sin(RAAN)));
    omega = wrapTo2Pi(u - nu);
end