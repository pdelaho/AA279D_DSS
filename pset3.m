%% Problem 1

close all; clc; clear;
mu = 398600.435436; % km^3/s^2

%% Part a: Chief and deputy orbital elements

% For the chief
a_c = 36943; % km, semi-major axis
e_c = 0.001; % eccentricity
inc_c = deg2rad(59); % inclination
omega_c = deg2rad(188); % argument of periapsis
RAAN_c = deg2rad(84); % right ascension of the ascending node
nu_c = deg2rad(0); % true anomaly
T = 2 * pi * sqrt(a_c^3 / mu);
n_c = sqrt(mu / a_c^3);

% For the deputy
delta_a = 0;
delta_e = - e_c * 0.001 * a_c;
delta_i = deg2rad(0.1) * a_c;
delta_omega = 0;
delta_RAAN = - deg2rad(0.03) * a_c;
delta_nu = deg2rad(-0.016) * a_c;

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

%% Functions

% function [pos_inertial, vel_inertial] = OEtoECI(a, e, inc, omega, RAAN, true_anom, mu)
%     cosE = (e + cos(true_anom)) / (1 + e * cos(true_anom));
%     sinE = sin(true_anom) * sqrt(1 - e^2) / (1 + e * cos(true_anom));
%     n = sqrt(mu / a^3);
%     
%     % Computing the position and velocity vectors in the perifocal frame
%     pos_perifocal = [a * (cosE - e) a * sqrt(1 - e^2) * sinE 0];
%     vel_perifocal = a * n / (1 - e * cosE) * [-sinE sqrt(1 - e^2)*cosE 0];
%     
%     % Computing the rotation matrix between perifocal and inertial frame
%     rotRAAN = [cos(RAAN) sin(-RAAN) 0;
%                -sin(-RAAN) cos(RAAN) 0;
%                0 0 1];
%     roti = [1 0 0;
%             0 cos(inc) sin(-inc);
%             0 -sin(-inc) cos(inc)];
%     rotomega = [cos(omega) sin(-omega) 0;
%                 -sin(-omega) cos(omega) 0;
%                 0 0 1];
%     rot_perifocalTOinertial = rotRAAN * roti * rotomega;
% 
%     % Rotating the position and velocity vectors
%     pos_inertial = rot_perifocalTOinertial * pos_perifocal';
%     vel_inertial = rot_perifocalTOinertial * vel_perifocal';
% end

function Rot = ECItoRTN(state)
    % Tring to enforce a shape
    % to have a stable output
    pos = zeros(1,3);
    vel = zeros(1,3);
    pos(:) = state(1:3);
    vel(:) = state(4:6);
    R = pos / norm(pos);
    h = cross(pos, vel);
    N = h / norm(h);
    T = cross(N, R);
    Rot = [R; T; N];
end