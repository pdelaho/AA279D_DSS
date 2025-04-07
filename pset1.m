%% Problem 2
close all; clc; clear;

path_config;
mu = 398600.435436; % km^3/s^2
R_E = 6378; % km
J2 = 0.108263e-2;

%% Part a: Orbital elements
% Change it to the actual values from the chosen mission
a = 36943; % km, semi-major axis
e = 0.8111; % eccentricity
inc = deg2rad(59); % inclination
omega = deg2rad(188); % argument of periapsis
RAAN = deg2rad(84); % right ascension of the ascending node
true_anom = 0; % true anomaly

%% Part b: Initial position and velocity in inertial frame
n = sqrt(mu / a^3); % mean motion
T = 2 * pi / n; % orbital period

[pos_inertial, vel_inertial] = OEtoECI(a, e, inc, omega, RAAN, true_anom, mu);

%% Part c: Orbit propagation with and without J2 effects
% Numerical propgation without J2 effects
init_state = [pos_inertial' vel_inertial'];
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 200);
N_unpert = 10000;
tspan = linspace(0, 10 * T, N_unpert);
[t_unpert, y_unpert] = ode89(@(t, state) state_der_posvel(t, state, mu), tspan, init_state, options);

% Numerical propagation with J2 effects
N_pert = 10000;
tspan = linspace(0, 10 * T, N_pert);
[t_pert, y_pert] = ode89(@(t, state) state_der_posvel_J2(t, state, mu, J2, R_E), tspan, init_state, options);

% Modeling the Earth
[X, Y, Z] = sphere(10);

% figure
% surf(X * R_E, Y * R_E, Z * R_E)
% hold on
% plot3(y_unpert(1:end, 1), y_unpert(1:end, 2), y_unpert(1:end, 3))
% plot3(y_pert(1:end, 1), y_pert(1:end, 2), y_pert(1:end, 3))
% axis equal
% grid on
% hold off
% legend('Earth', 'Unperturbed orbit', 'Perturbed orbit')
% xlabel('X-axis [km]')
% ylabel('Y-axis [km]')
% zlabel('Z-axis [km]')
% title('Unperturbed and J2-perturbed orbits around the Earth in ECI frame')
% 
% figure
% hold on
% plot(t_unpert / 3600, y_unpert(:, 1))
% plot(t_unpert / 3600, y_unpert(:, 2))
% plot(t_unpert / 3600, y_unpert(:, 3))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Position vector components [km]')
% legend('X-axis', 'Y-axis', 'Z-axis')
% title('Time evolution of the components of the position in ECI without J2 effects')
% 
% figure
% hold on
% plot(t_unpert / 3600, y_unpert(:, 4))
% plot(t_unpert / 3600, y_unpert(:, 5))
% plot(t_unpert / 3600, y_unpert(:, 6))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Velocity vector components [km/s]')
% legend('X-axis', 'Y-axis', 'Z-axis')
% title('Time evolution of the components of the velocity in ECI without J2 effects')
% 
% figure
% hold on
% plot(t_pert / 3600, y_pert(:, 1))
% plot(t_pert / 3600, y_pert(:, 2))
% plot(t_pert / 3600, y_pert(:, 3))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Position vector components [km]')
% legend('X-axis', 'Y-axis', 'Z-axis')
% title('Time evolution of the components of the position in ECI with J2 effects')
% 
% figure
% hold on
% plot(t_pert / 3600, y_pert(:, 4))
% plot(t_pert / 3600, y_pert(:, 5))
% plot(t_pert / 3600, y_pert(:, 6))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Velocity vector components [km/s]')
% legend('X-axis', 'Y-axis', 'Z-axis')
% title('Time evolution of the components of the velocity in ECI with J2 effects')

%% Part d: Analytical Keplerian propagation
M_unpert = n * t_unpert;
E_unpert = zeros(N_unpert, 1);
true_anom_unpert = zeros(N_unpert, 1);
y_unpert_Kep = zeros(size(y_unpert));
errors_RTN = zeros(size(y_unpert));

for i=1:N_unpert
    % Computing the eccentric anomaly
    E = eccentric_anom(M_unpert(i), e, 1e-10);
    E_unpert(i) = E;
    % Computing the true anomaly
    cos_nu = (cos(E) - e) / (1 - e * cos(E));
    sin_nu = sqrt(1 - e^2) * sin(E) / (1 - e * cos(E));
    nu = atan2(sin_nu, cos_nu);
    true_anom_unpert(i) = nu;
    % Computing the position and velocity vectors in the inertial frame
    [pos, vel] = OEtoECI(a, e, inc, omega, RAAN, nu, mu);
    y_unpert_Kep(i, :) = [pos' vel'];
    % Defining the RTN frame and rotation matrix
    R_vec = pos / norm(pos);
    N_vec = cross(pos, vel) / norm(cross(pos, vel));
    T_vec = cross(N_vec, R_vec);
    rot = [R_vec'; T_vec'; N_vec'];
    % Computing the angular velocity vector of inertial w.r.t. RTN
    theta_dot_vec = [0 0 -sqrt(mu / (a^3 * (1 - e^2)^3)) * (1 + e * cos(nu))^2];
    % Transforming the position and velocity vectors from inertial to RTN
    pos_RTN = rot * y_unpert(i, 1:3)';
    vel_RTN = rot * y_unpert(i, 4:6)' + cross(theta_dot_vec, pos_RTN)';
    pos_RTN_Kep = rot * y_unpert_Kep(i, 1:3)';
    vel_RTN_Kep = rot * y_unpert_Kep(i, 4:6)' + cross(theta_dot_vec, pos_RTN_Kep)';
    % Computing the errors in RTN frame
    error_pos = abs(pos_RTN - pos_RTN_Kep);
    error_vel = abs(vel_RTN - vel_RTN_Kep);
    errors_RTN(i, :) = [error_pos' error_vel'];
end

% figure
% hold on
% plot(t_unpert / 3600, errors_RTN(:, 1))
% plot(t_unpert / 3600, errors_RTN(:, 2))
% plot(t_unpert / 3600, errors_RTN(:, 3))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Error in absolute position [km]')
% title('Time evolution of the error in absolute position between numerical integration and analytical propagation in the RTN frame')
% legend('R-direction', 'T-direction', 'N-direction')
% 
% figure
% hold on
% plot(t_unpert / 3600, errors_RTN(:, 4))
% plot(t_unpert / 3600, errors_RTN(:, 5))
% plot(t_unpert / 3600, errors_RTN(:, 6))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Error in absolute velocity [km/s]')
% title('Time evolution of the error in absolute velocity between numerical integration and analytical propagation in the RTN frame')
% legend('R-direction', 'T-direction', 'N-direction')

%% Part e: Computing quantities throughout the orbit (un)perturbed
Kep_ele_unpert = zeros(N_unpert, 6);
ecc_vec_unpert = zeros(N_unpert, 3);
ang_mom_unpert = zeros(N_unpert, 3);
spec_energy_unpert = zeros(N_unpert, 1);

for i=1:N_unpert
    pos = y_unpert(i, 1:3);
    vel = y_unpert(i, 4:6);
    [a, e, inc, omega, RAAN, M] = Keplerian_elements(y_unpert(i,:), mu);
    Kep_ele_unpert(i, :) = [a, e, inc, omega, RAAN, M]';
    ecc_vec_unpert(i, :) = cross(vel, cross(pos, vel)) / mu - pos / norm(pos);
    % Or should we keep it in the perifocal coordinate?
    ang_mom_unpert(i, :) = cross(pos, vel);
    spec_energy_unpert(i) = norm(vel)^2 / 2 - mu / norm(pos);
end

% figure
% hold on
% plot(t_unpert / 3600, ecc_vec_unpert(:,1))
% plot(t_unpert / 3600, ecc_vec_unpert(:,2))
% plot(t_unpert / 3600, ecc_vec_unpert(:,3))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Components of the eccentricity vector [-]')
% legend('X-axis', 'Y-axis', 'Z-axis')
% title('Time variation of the eccentricity vector components in ECI for the unperturbed propagation')
% 
% figure
% hold on
% plot(t_unpert / 3600, ang_mom_unpert(:,1))
% plot(t_unpert / 3600, ang_mom_unpert(:,2))
% plot(t_unpert / 3600, ang_mom_unpert(:,3))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Components of the angular momentum vector [km^2/s]')
% legend('X-axis', 'Y-axis', 'Z-axis')
% title('Time variation of the angular momentum vector components in ECI for the unperturbed propagation')
% 
% figure
% plot(t_unpert / 3600, spec_energy_unpert)
% grid on
% xlabel('Time [hours]')
% ylabel('Specific mechanical energy [km^2/s^2]')
% title('Time variation of the specific mechanical energy for the unperturbed propagation')
% 
% figure
% plot(t_unpert / 3600, Kep_ele_unpert(:, 1))
% grid on
% xlabel('Time [hours]')
% ylabel('Semi-major axis [km]')
% title('Time variation of the semi-major axis for the unperturbed propagation')
% 
% figure
% plot(t_unpert / 3600, Kep_ele_unpert(:, 2))
% grid on
% xlabel('Time [hours]')
% ylabel('Eccentricity [-]')
% title('Time variation of the eccentricity for the unperturbed propagation')
% 
% % Maybe do different plots/subplots for each elements
% figure
% hold on
% plot(t_unpert / 3600, Kep_ele_unpert(:, 3:end))
% hold off
% grid on
% xlabel('Time [hours]')
% ylabel('Angular orbital elements [rad]')
% legend('Inclination', 'Argument of Periapsis', 'RAAN', 'Mean anomaly')
% title('Time variation of the inclination, argument of periapsis, RAAN and mean anomaly for the unperturbed propagation')

Kep_ele_pert = zeros(N_pert, 6);
ecc_vec_pert = zeros(N_pert, 3);
ang_mom_pert = zeros(N_pert, 3);
spec_energy_pert = zeros(N_pert, 1);

for i=1:N_pert
    pos = y_pert(i, 1:3);
    vel = y_pert(i, 4:6);
    [a, e, inc, omega, RAAN, M] = Keplerian_elements(y_pert(i,:), mu);
    Kep_ele_pert(i, :) = [a, e, inc, omega, RAAN, M]';
    ecc_vec_pert(i, :) = cross(vel, cross(pos, vel)) / mu - pos / norm(pos);
    ang_mom_pert(i, :) = cross(pos, vel);
    spec_energy_pert(i) = norm(vel)^2 / 2 - mu / norm(pos);
end

figure
hold on
plot(t_pert / 3600, ecc_vec_pert(:,1))
plot(t_pert / 3600, ecc_vec_pert(:,2))
plot(t_pert / 3600, ecc_vec_pert(:,3))
hold off
grid on
xlabel('Time [hours]')

title('Time variation of the eccentricity vector components for the perturbed propagation')

figure
hold on
plot(t_pert / 3600, ang_mom_pert(:,1))
plot(t_pert / 3600, ang_mom_pert(:,2))
plot(t_pert / 3600, ang_mom_pert(:,3))
hold off
grid on
xlabel('Time [hours]')

title('Time variation of the angular momentum vector components for the perturbed propagation')

figure
plot(t_pert / 3600, spec_energy_pert)
grid on
xlabel('Time [hours]')

title('Time variation of the specific mechanical energy for the perturbed propagation')

figure
plot(t_pert / 3600, Kep_ele_pert(:, 1))
grid on
xlabel('Time [hours]')

title('Time variation of the semi-major axis for the perturbed propagation')

figure
plot(t_pert / 3600, Kep_ele_pert(:, 2))
grid on
xlabel('Time [hours]')

title('Time variation of the eccentricity for the perturbed propagation')

figure
hold on
plot(t_pert / 3600, Kep_ele_pert(:, 3:end))
hold off
grid on
xlabel('Time [hours]')

legend('Inclination', 'Argument of periapsis', 'RAAN', 'Mean anomaly')
title('Time variation of the inclination, argument of periapsis, RAAN and mean anomaly for the perturbed propagation')

%% Part f: GVE with J2 effects
OE = [a * 1e3, e, inc, RAAN, omega, true_anom];
init_OE = osc2mean(OE);
a_mean = init_OE(1) / 1e3;
e_mean = init_OE(2);
inc_mean = init_OE(3);
RAAN_mean = init_OE(4);
omega_mean = init_OE(5);
M_mean = init_OE(6);

init_modified_OE = [a_mean e_mean * cos(omega_mean) e_mean * sin(omega_mean) inc_mean RAAN_mean omega_mean + true_anom];
tspan = linspace(0, 10 * T, N_pert);
[t_pert_OE, y_pert_OE] = ode89(@(t, state) GVEderJ2(t, state, mu, J2, R_E), tspan, init_modified_OE, options);

ecc_pert_GVE = zeros(N_pert, 1);
ecc_vec_pert_GVE = zeros(N_pert, 3);
omega_pert_GVE = zeros(N_pert, 1);
mean_anom_GVE = zeros(N_pert, 1);
ang_mom_pert_GVE = zeros(N_pert, 3);
spec_energy_pert_GVE = zeros(N_pert, 1);

for i=1:N_pert
    ecc = sqrt(y_pert_OE(i, 2)^2 + y_pert_OE(i, 3)^2);
    ecc_pert_GVE(i) = ecc;
    omega_pert = atan2(y_pert_OE(i, 3), y_pert_OE(i, 2));
    omega_pert_GVE(i) = omega_pert;
    true_anom_pert = y_pert_OE(i, 6) - omega_pert;
    cosE = (ecc + cos(true_anom_pert)) / (1 + ecc * cos(true_anom_pert));
    sinE = sin(true_anom_pert) * sqrt(1 - ecc^2) / (1 + ecc * cos(true_anom_pert));
    E = atan2(sinE, cosE);
    M = E - ecc * sinE;
    mean_anom_GVE(i) = M;
    [pos, vel] = OEtoECI(y_pert_OE(i, 1), ecc, y_pert_OE(i, 4), omega_pert, y_pert_OE(i, 5), true_anom_pert, mu);
    ang_mom_pert_GVE(i, :) = cross(pos, vel);
    spec_energy_pert_GVE(i) = - mu / (2 * y_pert_OE(i, 1));
    % eccentricity vector in ECI
    ecc_vec_pert_GVE(i, :) = cross(vel, cross(pos, vel)) / mu - pos / norm(pos);
end

% figure
% plot(t_pert_OE, y_pert_OE(:, 1))
% grid on
% title('Semi-major axis')
% 
% figure
% plot(t_pert_OE, y_pert_OE(:, 4))
% grid on
% title('Inclination')
% 
% figure
% plot(t_pert_OE, y_pert_OE(:, 5))
% grid on
% title('RAAN')
% 
% figure
% plot(t_pert_OE, omega_pert_GVE)
% grid on
% title('Argument of periapsis')
% 
% figure
% plot(t_pert_OE, ecc_pert_GVE)
% grid on
% title('Eccentricity')
% 
% figure
% plot(t_pert_OE, mean_anom_GVE)
% grid on
% title('Mean anomaly')
% 
% figure
% hold on
% plot(t_pert_OE, ang_mom_pert_GVE(:, 1))
% plot(t_pert_OE, ang_mom_pert_GVE(:, 2))
% plot(t_pert_OE, ang_mom_pert_GVE(:, 3))
% hold off
% title('Angular momentum')
% legend('x','y','z')
% 
% figure
% hold on
% plot(t_pert_OE, ecc_vec_pert_GVE(:, 1))
% plot(t_pert_OE, ecc_vec_pert_GVE(:, 2))
% plot(t_pert_OE, ecc_vec_pert_GVE(:, 3))
% hold off
% title('Eccentricity vector')
% legend('x', 'y', 'z')
% 
% figure
% plot(t_pert_OE, spec_energy_pert_GVE)
% grid on
% title('Specific mechanical energy')

%% Functions
function statedot = state_der_posvel(t, state, mu)
    statedot = zeros(size(state));
    statedot(1:3) = state(4:6);
    statedot(4:6) = - mu * state(1:3) / norm(state(1:3))^3;
end

function statedot = state_der_posvel_J2(t, state, mu, J2, R_E)
    statedot = zeros(size(state));
    statedot(1:3) = state(4:6);
    r = norm(state(1:3));
    d = - mu .* J2 .* R_E.^2 ./ 2 .* (6 .* state(3) ./ r.^5 .* [0 0 1] + (3 ./ r.^4 - 15 .* state(3).^2 ./ r.^6) .* state(1:3)' ./ r);
    statedot(4:6) = - mu .* state(1:3) ./ r.^3 + d';
end

function [a, e, i, omega, RAAN, M] = Keplerian_elements(state, mu)
    pos = state(1:3);
    vel = state(4:6);
    h = cross(pos, vel);
    W = h / norm(h);
    i = atan2(sqrt(W(1)^2 + W(2)^2), W(3)); % right atan?
    RAAN = atan2(W(1), - W(2));
    p = norm(h)^2 / mu;
    a = 1 / (2 / norm(pos) - norm(vel)^2 / mu);
    n = sqrt(mu / a^3);
    e = sqrt(1 - p / a);
    E = atan2((dot(pos, vel) / (a^2 * n)), (1 - norm(pos)/a));
    nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    M = E - e * sin(E);
    u = atan2((pos(3) / sin(i)), (pos(1) * cos(RAAN) + pos(2) * sin(RAAN)));
    omega = u - nu;
end

function E = eccentric_anom(M, e, epsilon)
    M = wrapTo2Pi(M);
    E = pi;
    while abs(- E + e * sin(E) + M) / (1 - e * cos(E)) > epsilon
        E_new = E - (E - e * sin(E) - M) / (1 - e * cos(E));
        E = E_new;
    end
end

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

function statedot = GVEderJ2(t, state, mu, J2, R_E)
    statedot = zeros(size(state));
    a = state(1);
    ex = state(2);
    ey = state(3);
    i = state(4);

    n = sqrt(mu / a^3);
    statedot(2) = - 3 / 4 * n * J2 * (R_E / (a * (1 - ex^2 - ey^2)))^2 * ey * (5 * cos(i)^2 - 1);
    statedot(3) = 3 / 4 * n * J2 * (R_E / (a * (1 - ex^2 - ey^2)))^2 * ex * (5 * cos(i)^2 - 1);
    statedot(5) = - 3 / 2 * n * J2 * (R_E / (a * (1 - ex^2 - ey^2)))^2 * cos(i);
    statedot(6) = 3 / 4 * n * J2 * (R_E / (a * (1 - ex^2 - ey^2)))^2 *(sqrt(1 - ex^2 - ey^2) * (3 * cos(i)^2 - 1) + 5 * cos(i)^2 - 1);
end