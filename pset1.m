%% Problem 2

close all; clc; clear;
path_config;

mu = 398600.435436; % km^3/s^2
R_E = 6378; % km
J2 = 0.108263e-2;

%% Part a: Orbital elements

a = 36943; % km
e = 0.8111;
inc = deg2rad(59);
omega = deg2rad(188);
RAAN = deg2rad(84);
f = 0;
n = sqrt(mu / a^3);
T = 2 * pi / n;

%% Part b: Initial position and velocity in inertial frame

[position_ECI, velocity_ECI] = OE2ECI(a, e, inc, omega, RAAN, f, mu);

%% Part c: Orbit propagation with and without J2 effects

% Numerical propgation without J2 effects
init_state = [position_ECI' velocity_ECI'];
options = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 100);
N = 10000;
tspan = linspace(0, 10 * T, N);
[t_unpert, y_unpert] = ode89(@(t, state) state_der_posvel(t, state, mu), tspan, init_state, options);

% Numerical propagation with J2 effects
[t_pert, y_pert] = ode89(@(t, state) state_der_posvel_J2(t, state, mu, J2, R_E), tspan, init_state, options);

% Modeling the Earth
[X, Y, Z] = sphere(10);

figure
surf(X * R_E, Y * R_E, Z * R_E)
hold on
plot3(y_unpert(1:end, 1), y_unpert(1:end, 2), y_unpert(1:end, 3))
plot3(y_pert(1:end, 1), y_pert(1:end, 2), y_pert(1:end, 3))
axis equal
grid on
hold off
legend('Earth', 'Unperturbed orbit', 'Perturbed orbit')
xlabel('X-axis [km]')
ylabel('Y-axis [km]')
zlabel('Z-axis [km]')

figure
subplot(3,2,1)
hold on
plot(t_unpert / T, y_unpert(:, 1))
plot(t_pert / T, y_pert(:, 1))
hold off
grid on
legend('No J2', 'J2')
ylabel('X-axis component [km]')
title('Position')

subplot(3,2,3)
hold on
plot(t_unpert / T, y_unpert(:, 2))
plot(t_pert / T, y_pert(:, 2))
hold off
grid on
ylabel('Y-axis component [km]')

subplot(3,2,5)
hold on
plot(t_unpert / T, y_unpert(:, 3))
plot(t_pert / T, y_pert(:, 3))
hold off
grid on
ylabel('Z-axis component [km]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_unpert / T, y_unpert(:, 4))
plot(t_pert / T, y_pert(:, 4))
hold off
grid on
ylabel('X-axis component [km/s]')
title('Velocity')

subplot(3,2,4)
hold on
plot(t_unpert / T, y_unpert(:, 5))
plot(t_pert / T, y_pert(:, 5))
hold off
grid on
ylabel('Y-axis component [km/s]')

subplot(3,2,6)
hold on
plot(t_unpert / T, y_unpert(:, 6))
plot(t_pert / T, y_pert(:, 6))
hold off
grid on
ylabel('Z-axis component [km/s]')
xlabel('Orbital Periods')

%% Part d: Analytical Keplerian propagation

M_unpert = n * t_unpert;
errors_RTN = zeros(size(y_unpert));

for i=1:N
    % Computing the eccentric and true anomalies
    E = M2E(M_unpert(i), e, 1e-12);
    f = E2f(E, e);

    % Computing the position and velocity vectors in the inertial frame
    [pos, vel] = OE2ECI(a, e, inc, omega, RAAN, f, mu);

    % Computing the errors in ECI frame
    error_pos_ECI = y_unpert(i, 1:3)' - pos;
    error_vel_ECI = y_unpert(i, 4:6)' - vel;
    
    % Transforming the errors from ECI to RTN
    state_RTN = ECI2RTN([pos', vel'], [error_pos_ECI', error_vel_ECI'], mu);
    errors_RTN(i, :) = abs(state_RTN);
end

figure
subplot(3,2,1)
plot(t_unpert / T, errors_RTN(:, 1))
grid on
ylabel('R-direction [km]')
title('Error in position')

subplot(3,2,3)
plot(t_unpert / T, errors_RTN(:, 2))
grid on
ylabel('T-direction [km]')

subplot(3,2,5)
plot(t_unpert / T, errors_RTN(:, 3))
grid on
ylabel('N-direction [km]')
xlabel('Orbital Periods')

subplot(3,2,2)
plot(t_unpert / T, errors_RTN(:, 4))
grid on
ylabel('R-direction [km/s]')
title('Error in velocity')

subplot(3,2,4)
plot(t_unpert / T, errors_RTN(:, 5))
grid on
ylabel('T-direction [km/s]')

subplot(3,2,6)
plot(t_unpert / T, errors_RTN(:, 6))
grid on
ylabel('N-direction [km/s]')
xlabel('Orbital Periods')

%% Part e: Computing quantities throughout the orbit (un)perturbed

OE_unpert = zeros(N, 6);
evec_unpert = zeros(N, 3);
ang_mom_unpert = zeros(N, 3);
spec_energy_unpert = zeros(N, 1);

OE_pert = zeros(N, 6);
evec_pert = zeros(N, 3);
ang_mom_pert = zeros(N, 3);
spec_energy_pert = zeros(N, 1);

for i=1:N
    % No J2 perturbations
    pos = y_unpert(i, 1:3);
    vel = y_unpert(i, 4:6);
    [a_unpert, e_unpert, inc_unpert, omega_unpert, RAAN_unpert, M_unpert] = ECI2OE_M(y_unpert(i, :), mu);
    OE_unpert(i, :) = [a_unpert, e_unpert, inc_unpert, omega_unpert, RAAN_unpert, M_unpert];
    evec_unpert(i, :) = cross(vel, cross(pos, vel)) / mu - pos / norm(pos);
    ang_mom_unpert(i, :) = cross(pos, vel);
    spec_energy_unpert(i) = norm(vel)^2 / 2 - mu / norm(pos);

    % J2 perturbations
    pos = y_pert(i, 1:3);
    vel = y_pert(i, 4:6);
    [a_pert, e_pert, inc_pert, omega_pert, RAAN_pert, M_pert] = ECI2OE_M(y_pert(i,:), mu);
    OE_pert(i, :) = [a_pert, e_pert, inc_pert, omega_pert, RAAN_pert, M_pert]';
    evec_pert(i, :) = cross(vel, cross(pos, vel)) / mu - pos / norm(pos);
    ang_mom_pert(i, :) = cross(pos, vel);
    spec_energy_pert(i) = norm(vel)^2 / 2 - mu / norm(pos);
end

figure
subplot(3,2,1)
hold on
plot(t_unpert / T, OE_unpert(:, 1))
plot(t_pert / T, OE_pert(:, 1))
hold off
grid on
legend('No J2', 'J2')
ylabel('Semi-major axis [km]')

subplot(3,2,2)
hold on
plot(t_unpert / T, OE_unpert(:, 4))
plot(t_pert / T, OE_pert(:, 4))
hold off
grid on
ylabel('Argument of periapsis [rad]')

subplot(3,2,3)
hold on
plot(t_unpert / T, OE_unpert(:, 2))
plot(t_pert / T, OE_pert(:, 2))
hold off
grid on
ylabel('Eccentricity [-]')

subplot(3,2,4)
hold on
plot(t_unpert / T, OE_unpert(:, 5))
plot(t_pert / T, OE_pert(:, 5))
hold off
grid on
ylabel('RAAN [rad]')

subplot(3,2,5)
hold on
plot(t_unpert / T, OE_unpert(:, 3))
plot(t_pert / T, OE_pert(:, 3))
hold off
grid on
ylabel('Inclination [rad]')
xlabel('Orbital Periods')

subplot(3,2,6)
hold on
plot(t_unpert / T, OE_unpert(:, 6))
plot(t_pert / T, OE_pert(:, 6))
hold off
grid on
ylabel('Mean anomaly [rad]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(t_unpert / T, evec_unpert(:,1))
plot(t_pert / T, evec_pert(:,1))
hold off
grid on
legend('No J2', 'J2')
ylabel('X-component [-]')

subplot(3,1,2)
hold on
plot(t_unpert / T, evec_unpert(:,2))
plot(t_pert / T, evec_pert(:,2))
hold off
grid on
ylabel('Y-component [-]')

subplot(3,1,3)
hold on
plot(t_unpert / T, evec_unpert(:,3))
plot(t_pert / T, evec_pert(:,3))
hold off
grid on
ylabel('Z-component [-]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(t_unpert / T, ang_mom_unpert(:,1))
plot(t_pert / T, ang_mom_pert(:,1))
hold off
grid on
legend('No J2', 'J2')
ylabel('X-component [km^2/s]')

subplot(3,1,2)
hold on
plot(t_unpert / T, ang_mom_unpert(:,2))
plot(t_pert / T, ang_mom_pert(:,2))
hold off
grid on
ylabel('Y-component [km^2/s]')

subplot(3,1,3)
hold on
plot(t_unpert / T, ang_mom_unpert(:,3))
plot(t_pert / T, ang_mom_pert(:,3))
hold off
grid on
ylabel('Z-component [km^2/s]')
xlabel('Orbital Periods')

figure
hold on
plot(t_unpert / T, spec_energy_unpert)
plot(t_pert / T, spec_energy_pert)
hold off
grid on
ylabel('Specific mechanical energy [km^2/s^2]')
xlabel('Orbital Periods')
legend('No J2', 'J2')

%% Part f: GVE with J2 effects

initial_osculating_OE = [a e * cos(omega) e * sin(omega) inc RAAN omega + 0];
[t_osculating_OE, y_osculating_OE] = ode89(@(t, state) GVEderJ2(t, state, mu, J2, R_E), tspan, initial_osculating_OE, options);

e_osculating_GVE = zeros(N, 1);
evec_osculating_GVE = zeros(N, 3);
omega_osculating_GVE = zeros(N, 1);
M_osculating_GVE = zeros(N, 1);
ang_mom_osculating_GVE = zeros(N, 3);
spec_energy_osculating_GVE = zeros(N, 1);

for i=1:N
    ecc = sqrt(y_osculating_OE(i, 2)^2 + y_osculating_OE(i, 3)^2);
    e_osculating_GVE(i) = ecc;
    omega_pert = atan2(y_osculating_OE(i, 3), y_osculating_OE(i, 2));
    omega_osculating_GVE(i) = omega_pert;
    M_pert = y_osculating_OE(i, 6) - omega_pert;
    M_osculating_GVE(i) = M_pert;
    E_pert = M2E(M_pert, ecc, 1e-12);
    f_pert = E2f(E_pert, ecc);
    [pos, vel] = OE2ECI(y_osculating_OE(i, 1), ecc, y_osculating_OE(i, 4), omega_pert, y_osculating_OE(i, 5), f_pert, mu);
    ang_mom_osculating_GVE(i, :) = cross(pos, vel);
    spec_energy_osculating_GVE(i) = - mu / (2 * y_osculating_OE(i, 1));
    evec_osculating_GVE(i, :) = cross(vel, cross(pos, vel)) / mu - pos / norm(pos);
end

figure
subplot(3,2,1)
hold on
plot(t_osculating_OE / T, OE_pert(:, 1))
plot(t_osculating_OE / T, y_osculating_OE(:, 1))
hold off
grid on
legend('Osculating', 'Mean')
ylabel('Semi-majo axis [km]')

subplot(3,2,3)
hold on
plot(t_osculating_OE / T, OE_pert(:, 2))
plot(t_osculating_OE / T, e_osculating_GVE)
hold off
grid on
ylabel('Eccentricity [-]')

subplot(3,2,5)
hold on
plot(t_osculating_OE / T, OE_pert(:, 3))
plot(t_osculating_OE / T, y_osculating_OE(:, 4))
hold off
grid on
ylabel('Inclination [rad]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_osculating_OE / T, wrapToPi(OE_pert(:, 4)))
plot(t_osculating_OE / T, wrapToPi(omega_osculating_GVE))
hold off
grid on
ylabel('Argument of periapsis [rad]')

subplot(3,2,4)
hold on
plot(t_osculating_OE / T, OE_pert(:, 5))
plot(t_osculating_OE / T, y_osculating_OE(:, 5))
hold off
grid on
ylabel('RAAN [rad]')

subplot(3,2,6)
hold on
plot(t_osculating_OE / T, OE_pert(:, 6))
plot(t_osculating_OE / T, wrapTo2Pi(M_osculating_GVE))
hold off
grid on
ylabel('Mean anomaly [rad]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(t_osculating_OE / T, evec_pert(:, 1))
plot(t_osculating_OE / T, evec_osculating_GVE(:, 1))
hold off
grid on
legend('Osculating', 'Mean')
ylabel('X-component [-]')

subplot(3,1,2)
hold on
plot(t_osculating_OE / T, evec_pert(:, 2))
plot(t_osculating_OE / T, evec_osculating_GVE(:, 2))
hold off
grid on
ylabel('Y-component [-]')

subplot(3,1,3)
hold on
plot(t_osculating_OE / T, evec_pert(:, 3))
plot(t_osculating_OE / T, evec_osculating_GVE(:, 3))
hold off
grid on
ylabel('Z-component [-]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(t_osculating_OE / T, ang_mom_pert(:, 1))
plot(t_osculating_OE / T, ang_mom_osculating_GVE(:, 1))
hold off
grid on
legend('Osculating', 'Mean')
ylabel('X-component [km^2/s]')

subplot(3,1,2)
hold on
plot(t_osculating_OE / T, ang_mom_pert(:, 2))
plot(t_osculating_OE / T, ang_mom_osculating_GVE(:, 2))
hold off
grid on
ylabel('Y-component [km^2/s]')

subplot(3,1,3)
hold on
plot(t_osculating_OE / T, ang_mom_pert(:, 3))
plot(t_osculating_OE / T, ang_mom_osculating_GVE(:, 3))
hold off
grid on
ylabel('Z-component [km^2/s]')
xlabel('Orbital Periods')

figure
hold on
plot(t_osculating_OE / T, spec_energy_pert)
plot(t_osculating_OE / T, spec_energy_osculating_GVE)
hold off
grid on
legend('Osculating', 'Mean')
ylabel('Specific mechanical energy [km^2/s^2]')
xlabel('Orbital Periods')

%% Part g: GVE with mean initial elements

OE = [a * 1e3, e, inc, RAAN, omega, f];
init_OE = osc2mean(OE);
a_mean = init_OE(1) / 1e3;
e_mean = init_OE(2);
inc_mean = init_OE(3);
RAAN_mean = init_OE(4);
omega_mean = init_OE(5);
M_mean = init_OE(6);

init_modified_OE = [a_mean e_mean * cos(omega_mean) e_mean * sin(omega_mean) inc_mean RAAN_mean omega_mean + M_mean];
[t_mean_OE, y_mean_OE] = ode89(@(t, state) GVEderJ2(t, state, mu, J2, R_E), tspan, init_modified_OE, options);

e_mean_GVE = zeros(N, 1);
evec_mean_GVE = zeros(N, 3);
omega_pert_GVE = zeros(N, 1);
M_mean_GVE = zeros(N, 1);
ang_mom_mean_GVE = zeros(N, 3);
spec_energy_mean_GVE = zeros(N, 1);

for i=1:N
    ecc = sqrt(y_mean_OE(i, 2)^2 + y_mean_OE(i, 3)^2);
    e_mean_GVE(i) = ecc;
    omega_pert = atan2(y_mean_OE(i, 3), y_mean_OE(i, 2));
    omega_pert_GVE(i) = omega_pert;
    M_pert = y_mean_OE(i, 6) - omega_pert;
    M_mean_GVE(i) = M_pert;
    E_pert = M2E(M_pert, ecc, 1e-12);
    f_pert = E2f(E_pert, ecc);
    [pos, vel] = OE2ECI(y_mean_OE(i, 1), ecc, y_mean_OE(i, 4), omega_pert, y_mean_OE(i, 5), f_pert, mu);
    ang_mom_mean_GVE(i, :) = cross(pos, vel);
    spec_energy_mean_GVE(i) = - mu / (2 * y_mean_OE(i, 1));
    evec_mean_GVE(i, :) = cross(vel, cross(pos, vel)) / mu - pos / norm(pos);
end

figure
subplot(3,2,1)
hold on
plot(t_mean_OE / T, OE_pert(:, 1))
plot(t_mean_OE / T, y_mean_OE(:, 1))
hold off
grid on
legend('Osculating', 'Mean')
ylabel('Semi-majo axis [km]')

subplot(3,2,3)
hold on
plot(t_mean_OE / T, OE_pert(:, 2))
plot(t_mean_OE / T, e_mean_GVE)
hold off
grid on
ylabel('Eccentricity [-]')

subplot(3,2,5)
hold on
plot(t_mean_OE / T, OE_pert(:, 3))
plot(t_mean_OE / T, y_mean_OE(:, 4))
hold off
grid on
ylabel('Inclination [rad]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_mean_OE / T, wrapToPi(OE_pert(:, 4)))
plot(t_mean_OE / T, wrapToPi(omega_pert_GVE))
hold off
grid on
ylabel('Argument of periapsis [rad]')

subplot(3,2,4)
hold on
plot(t_mean_OE / T, OE_pert(:, 5))
plot(t_mean_OE / T, y_mean_OE(:, 5))
hold off
grid on
ylabel('RAAN [rad]')

subplot(3,2,6)
hold on
plot(t_mean_OE / T, OE_pert(:, 6))
plot(t_mean_OE / T, wrapTo2Pi(M_mean_GVE))
hold off
grid on
ylabel('Mean anomaly [rad]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(t_mean_OE / T, evec_pert(:, 1))
plot(t_mean_OE / T, evec_mean_GVE(:, 1))
hold off
grid on
legend('Osculating', 'Mean')
ylabel('X-component [-]')

subplot(3,1,2)
hold on
plot(t_mean_OE / T, evec_pert(:, 2))
plot(t_mean_OE / T, evec_mean_GVE(:, 2))
hold off
grid on
ylabel('Y-component [-]')

subplot(3,1,3)
hold on
plot(t_mean_OE / T, evec_pert(:, 3))
plot(t_mean_OE / T, evec_mean_GVE(:, 3))
hold off
grid on
ylabel('Z-component [-]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(t_mean_OE / T, ang_mom_pert(:, 1))
plot(t_mean_OE / T, ang_mom_mean_GVE(:, 1))
hold off
grid on
legend('Osculating', 'Mean')
ylabel('X-component [km^2/s]')

subplot(3,1,2)
hold on
plot(t_mean_OE / T, ang_mom_pert(:, 2))
plot(t_mean_OE / T, ang_mom_mean_GVE(:, 2))
hold off
grid on
ylabel('Y-component [km^2/s]')

subplot(3,1,3)
hold on
plot(t_mean_OE / T, ang_mom_pert(:, 3))
plot(t_mean_OE / T, ang_mom_mean_GVE(:, 3))
hold off
grid on
ylabel('Z-component [km^2/s]')
xlabel('Orbital Periods')

figure
hold on
plot(t_mean_OE / T, spec_energy_pert)
plot(t_mean_OE / T, spec_energy_mean_GVE)
hold off
grid on
legend('Osculating', 'Mean')
ylabel('Specific mechanical energy [km^2/s^2]')
xlabel('Orbital Periods')

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
    statedot(6) = 3 / 4 * n * J2 * (R_E / (a * (1 - ex^2 - ey^2)))^2 *(sqrt(1 - ex^2 - ey^2) * (3 * cos(i)^2 - 1) + 5 * cos(i)^2 - 1) + n;
end