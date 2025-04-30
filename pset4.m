%% Problem 1

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
R_E = 6378; % km
J2 = 0.108263e-2;

%% Question 1: Chief's initial condition

a_chief = 36943; % km
e_chief = 0.0001;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
nu_chief = 0;
M_chief = 0;
T = 2 * pi * sqrt(a_chief^3 / mu);
n_chief = sqrt(mu / a_chief^3);
init_oe_chief = [a_chief * 1e3, e_chief, inc_chief, RAAN_chief, omega_chief, M_chief];
init_mean_oe_chief = osc2mean(init_oe_chief);
initial_mean_oe_chief = [init_mean_oe_chief(1:3)', init_mean_oe_chief(5), init_mean_oe_chief(4), init_mean_oe_chief(6)];

%% Question 2: Deputy's initial condition

delta_a = 0;
delta_lambda = 100e-3 / a_chief;
delta_ex = 50e-3 / a_chief;
delta_ey = 100e-3 / a_chief;
delta_ix = 30e-3 / a_chief;
delta_iy = 200e-3 / a_chief;

a_deputy = a_chief + delta_a * a_chief;
inc_deputy = inc_chief + delta_ix;
RAAN_deputy = RAAN_chief + delta_iy / sin(inc_chief);
e_deputy = sqrt((delta_ex + e_chief * cos(omega_chief))^2 + (delta_ey + e_chief * sin(omega_chief))^2);
omega_deputy = atan2(delta_ey + e_chief * sin(omega_chief), delta_ex + e_chief * cos(omega_chief));
M_deputy = delta_lambda + M_chief + omega_chief - (RAAN_deputy - RAAN_chief) * cos(inc_chief) - omega_deputy;
E_deputy = M2E(M_deputy, e_deputy, 1e-12);
nu_deputy = E2f(E_deputy, e_deputy);
n_deputy = sqrt(mu / a_deputy^3);
init_oe_deputy = [a_deputy * 1e3, e_deputy, inc_deputy, RAAN_deputy, omega_deputy, M_deputy];
init_mean_oe_deputy = osc2mean(init_oe_deputy);
initial_mean_oe_deputy = [init_mean_oe_deputy(1:3)', init_mean_oe_deputy(5), init_mean_oe_deputy(4), init_mean_oe_deputy(6)];

%% Question 3: Numerical integration of FODE

[pos_chief, vel_chief] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief, mu);
[pos_deputy_1, vel_deputy_1] = OE2ECI(a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, nu_deputy, mu);

init_cond_1 = [pos_chief', vel_chief', pos_deputy_1', vel_deputy_1'];

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100);
N = 15000;
tspan = linspace(0, 15 * T, N);

% Numerical integration without J2 perturbations
[t_unpert_1, y_unpert_1] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, init_cond_1, options);

oe_c_osc_unpert = zeros(N, 6);
oe_d_osc_unpert = zeros(N, 6);
roe_osc_unpert = zeros(N, 6);

oe_c_mean_unpert = zeros(N, 6);
oe_d_mean_unpert = zeros(N, 6);
roe_mean_unpert = zeros(N, 6);

for i=1:N
    [a_c, e_c, i_c, omega_c, RAAN_c, M_c] = ECI2OE_M(y_unpert_1(i,1:6), mu);
    oe_c = [a_c, e_c, i_c, omega_c, RAAN_c, M_c];
    [a_d, e_d, i_d, omega_d, RAAN_d, M_d] = ECI2OE_M(y_unpert_1(i,7:12), mu);
    oe_d = [a_d, e_d, i_d, omega_d, RAAN_d, M_d];
    qnsoe_c = OE2QNSOE(oe_c);
    qnsoe_d = OE2QNSOE(oe_d);
    oe_c_osc_unpert(i, :) = qnsoe_c;
    oe_d_osc_unpert(i, :) = qnsoe_d;
    ROE = OE2ROE(oe_c, oe_d);
    roe_osc_unpert(i, :) = ROE;

    oe_c_mean = osc2mean([oe_c(1)*1e3, oe_c(2:3), oe_c(5), oe_c(4), oe_c(6)], 0);
    oe_d_mean = osc2mean([oe_d(1)*1e3, oe_d(2:3), oe_d(5), oe_d(4), oe_d(6)], 0);
    oe_c_mean = [oe_c_mean(1)*1e-3, oe_c_mean(2:3)', oe_c_mean(5), oe_c_mean(4), oe_c_mean(6)];
    oe_d_mean = [oe_d_mean(1)*1e-3, oe_d_mean(2:3)', oe_d_mean(5), oe_d_mean(4), oe_d_mean(6)];
    qnsoe_c_mean = OE2QNSOE(oe_c_mean);
    qnsoe_d_mean = OE2QNSOE(oe_d_mean);
    oe_c_mean_unpert(i, :) = qnsoe_c_mean;
    oe_d_mean_unpert(i, :) = qnsoe_d_mean;
    ROE_mean = OE2ROE(oe_c_mean, oe_d_mean);
    roe_mean_unpert(i, :) = ROE_mean;
end

figure
subplot(3,2,1)
hold on
plot(t_unpert_1 / T, oe_c_osc_unpert(:, 1), 'b-', 'LineWidth', 3)
plot(t_unpert_1 / T, oe_c_mean_unpert(:, 1), 'r--', 'LineWidth', 3)
hold off
grid on
legend('Osculating', 'Mean')
ylabel('Semi-major axis [km]')

subplot(3,2,3)
hold on
plot(t_unpert_1 / T, oe_c_osc_unpert(:, 2), 'b-', 'LineWidth', 2)
plot(t_unpert_1 / T, oe_c_mean_unpert(:, 2), 'r--', 'LineWidth', 2)
hold off
grid on
ylabel('Mean argument of latitude [rad]')

subplot(3,2,5)
hold on
plot(t_unpert_1 / T, oe_c_osc_unpert(:, 3), 'b-', 'LineWidth', 2)
plot(t_unpert_1 / T, oe_c_mean_unpert(:, 3), 'r--', 'LineWidth', 2)
hold off
grid on
ylabel('e_x [-]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_unpert_1 / T, oe_c_osc_unpert(:, 4), 'b-', 'LineWidth', 2)
plot(t_unpert_1 / T, oe_c_mean_unpert(:, 4), 'r--', 'LineWidth', 2)
hold off
grid on
ylabel('e_y [-]')

subplot(3,2,4)
hold on
plot(t_unpert_1 / T, oe_c_osc_unpert(:, 5), 'b-', 'LineWidth', 3)
plot(t_unpert_1 / T, oe_c_mean_unpert(:, 5), 'r--', 'LineWidth', 3)
hold off
grid on
ylabel('Inclination [rad]')

subplot(3,2,6)
hold on
plot(t_unpert_1 / T, oe_c_osc_unpert(:, 6), 'b-', 'LineWidth', 3)
plot(t_unpert_1 / T, oe_c_mean_unpert(:, 6), 'r--', 'LineWidth', 3)
hold off
grid on
ylabel('RAAN [rad]')
xlabel('Orbital Periods')

figure
subplot(3,2,1)
hold on
plot(t_unpert_1 / T, oe_d_osc_unpert(:, 1), 'b-', 'LineWidth', 3)
plot(t_unpert_1 / T, oe_d_mean_unpert(:, 1), 'r--', 'LineWidth', 3)
hold off
grid on
legend('Osculating', 'Mean')
ylabel('Semi-major axis [km]')

subplot(3,2,3)
hold on
plot(t_unpert_1 / T, oe_d_osc_unpert(:, 2), 'b-', 'LineWidth', 2)
plot(t_unpert_1 / T, oe_d_mean_unpert(:, 2), 'r--', 'LineWidth', 2)
hold off
grid on
ylabel('Mean argument of latitude [rad]')

subplot(3,2,5)
hold on
plot(t_unpert_1 / T, oe_d_osc_unpert(:, 3), 'b-', 'LineWidth', 2)
plot(t_unpert_1 / T, oe_d_mean_unpert(:, 3), 'r--', 'LineWidth', 2)
hold off
grid on
ylabel('e_x [-]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_unpert_1 / T, oe_d_osc_unpert(:, 4), 'b-', 'LineWidth', 2)
plot(t_unpert_1 / T, oe_d_mean_unpert(:, 4), 'r--', 'LineWidth', 2)
hold off
grid on
ylabel('e_y [-]')

subplot(3,2,4)
hold on
plot(t_unpert_1 / T, oe_d_osc_unpert(:, 5), 'b-', 'LineWidth', 3)
plot(t_unpert_1 / T, oe_d_mean_unpert(:, 5), 'r--', 'LineWidth', 3)
hold off
grid on
ylabel('Inclination [rad]')

subplot(3,2,6)
hold on
plot(t_unpert_1 / T, oe_d_osc_unpert(:, 6), 'b-', 'LineWidth', 3)
plot(t_unpert_1 / T, oe_d_mean_unpert(:, 6), 'r--', 'LineWidth', 3)
hold off
grid on
ylabel('RAAN [rad]')
xlabel('Orbital Periods')

figure
subplot(3,2,1)
hold on
plot(t_unpert_1 / T, roe_osc_unpert(:, 1) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert(:, 1) * a_chief * 1e3, 'r--')
hold off
grid on
legend('Osculating', 'Mean')
ylabel('a_c \deltaa [m]')

subplot(3,2,3)
hold on
plot(t_unpert_1 / T, roe_osc_unpert(:, 2) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert(:, 2) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('a_c \delta\lambda [m]')

subplot(3,2,5)
hold on
plot(t_unpert_1 / T, roe_osc_unpert(:, 3) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert(:, 3) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('a_c \deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_unpert_1 / T, roe_osc_unpert(:, 4) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert(:, 4) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('a_c \deltae_y [m]')

subplot(3,2,4)
hold on
plot(t_unpert_1 / T, roe_osc_unpert(:, 5) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert(:, 5) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('a_c \deltai_x [m]')

subplot(3,2,6)
hold on
plot(t_unpert_1 / T, roe_osc_unpert(:, 6) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert(:, 6) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('a_c \deltai_y [m]')
xlabel('Orbital Periods')

% Numerical integration with J2 perturbations
[t_pert_1, y_pert_1] = ode89(@(t, state) state_der_posvel_J2_2sats(t, state, mu, J2, R_E), tspan, init_cond_1, options);

oe_c_osc_pert = zeros(N, 6);
oe_d_osc_pert = zeros(N, 6);
roe_osc_pert = zeros(N, 6);

oe_c_mean_pert = zeros(N, 6);
oe_d_mean_pert = zeros(N, 6);
roe_mean_pert = zeros(N, 6);

u = zeros(N,1);

for i=1:N
    [a_c, e_c, i_c, omega_c, RAAN_c, M_c] = ECI2OE_M(y_pert_1(i,1:6), mu);
%     [a_c, e_c, i_c, omega_c, RAAN_c, M_c] = inert2OE(y_pert_1(i,1:6), mu);
    oe_c = [a_c, e_c, i_c, omega_c, RAAN_c, M_c];
    [a_d, e_d, i_d, omega_d, RAAN_d, M_d] = ECI2OE_M(y_pert_1(i,7:12), mu);
%     [a_d, e_d, i_d, omega_d, RAAN_d, M_d] = inert2OE(y_pert_1(i,7:12), mu);
    oe_d = [a_d, e_d, i_d, omega_d, RAAN_d, M_d];
    u(i) = M_c+omega_c;
    qnsoe_c = OE2QNSOE(oe_c);
    qnsoe_d = OE2QNSOE(oe_d);
    oe_c_osc_pert(i, :) = qnsoe_c;
    oe_d_osc_pert(i, :) = qnsoe_d;
    ROE = OE2ROE(oe_c, oe_d);
    roe_osc_pert(i, :) = ROE;

    oe_c_mean = osc2mean([oe_c(1)*1e3, oe_c(2:3), oe_c(5), oe_c(4), oe_c(6)]);
    oe_d_mean = osc2mean([oe_d(1)*1e3, oe_d(2:3), oe_d(5), oe_d(4), oe_d(6)]);
    oe_c_mean = [oe_c_mean(1)*1e-3, oe_c_mean(2:3)', oe_c_mean(5), oe_c_mean(4), oe_c_mean(6)];
    oe_d_mean = [oe_d_mean(1)*1e-3, oe_d_mean(2:3)', oe_d_mean(5), oe_d_mean(4), oe_d_mean(6)];
    qnsoe_c_mean = OE2QNSOE(oe_c_mean);
    qnsoe_d_mean = OE2QNSOE(oe_d_mean);
    oe_c_mean_pert(i, :) = qnsoe_c_mean;
    oe_d_mean_pert(i, :) = qnsoe_d_mean;
    ROE_mean = OE2ROE(oe_c_mean, oe_d_mean);
    roe_mean_pert(i, :) = ROE_mean;
end

figure
subplot(3,2,1)
hold on
plot(t_pert_1 / T, oe_c_osc_pert(:, 1))
plot(t_pert_1 / T, oe_c_mean_pert(:, 1))
hold off
grid on
legend('Osculating', 'Mean')
ylabel('Semi-major axis [km]')

subplot(3,2,3)
hold on
plot(t_pert_1 / T, oe_c_osc_pert(:, 2))
plot(t_pert_1 / T, oe_c_mean_pert(:, 2))
hold off
grid on
ylabel('Mean argument of latitude [rad]')

subplot(3,2,5)
hold on
plot(t_pert_1 / T, oe_c_osc_pert(:, 3))
plot(t_pert_1 / T, oe_c_mean_pert(:, 3))
hold off
grid on
ylabel('e_x [-]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_pert_1 / T, oe_c_osc_pert(:, 4))
plot(t_pert_1 / T, oe_c_mean_pert(:, 4))
hold off
grid on
ylabel('e_y [-]')

subplot(3,2,4)
hold on
plot(t_pert_1 / T, oe_c_osc_pert(:, 5))
plot(t_pert_1 / T, oe_c_mean_pert(:, 5))
hold off
grid on
ylabel('Inclination [rad]')

subplot(3,2,6)
hold on
plot(t_pert_1 / T, oe_c_osc_pert(:, 6))
plot(t_pert_1 / T, oe_c_mean_pert(:, 6))
hold off
grid on
ylabel('RAAN [rad]')
xlabel('Orbital Periods')

figure
subplot(3,2,1)
hold on
plot(t_pert_1 / T, oe_d_osc_pert(:, 1))
plot(t_pert_1 / T, oe_d_mean_pert(:, 1))
hold off
grid on
legend('Osculating', 'Mean')
ylabel('Semi-major axis [km]')

subplot(3,2,3)
hold on
plot(t_pert_1 / T, oe_d_osc_pert(:, 2))
plot(t_pert_1 / T, oe_d_mean_pert(:, 2))
hold off
grid on
ylabel('Mean argument of latitude [rad]')

subplot(3,2,5)
hold on
plot(t_pert_1 / T, oe_d_osc_pert(:, 3))
plot(t_pert_1 / T, oe_d_mean_pert(:, 3))
hold off
grid on
ylabel('e_x [-]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_pert_1 / T, oe_d_osc_pert(:, 4))
plot(t_pert_1 / T, oe_d_mean_pert(:, 4))
hold off
grid on
ylabel('e_y [-]')

subplot(3,2,4)
hold on
plot(t_pert_1 / T, oe_d_osc_pert(:, 5))
plot(t_pert_1 / T, oe_d_mean_pert(:, 5))
hold off
grid on
ylabel('Inclination [rad]')

subplot(3,2,6)
hold on
plot(t_pert_1 / T, oe_d_osc_pert(:, 6))
plot(t_pert_1 / T, oe_d_mean_pert(:, 6))
hold off
grid on
ylabel('RAAN [rad]')
xlabel('Orbital Periods')

figure
subplot(3,2,1)
hold on
plot(t_pert_1 / T, roe_osc_pert(:, 1) * a_chief * 1e3)
plot(t_pert_1 / T, roe_mean_pert(:, 1) * a_chief * 1e3)
hold off
grid on
legend('Osculating', 'Mean')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(t_pert_1 / T, roe_osc_pert(:, 2) * a_chief * 1e3)
plot(t_pert_1 / T, roe_mean_pert(:, 2) * a_chief * 1e3)
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(t_pert_1 / T, roe_osc_pert(:, 3) * a_chief * 1e3)
plot(t_pert_1 / T, roe_mean_pert(:, 3) * a_chief * 1e3)
hold off
grid on
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_pert_1 / T, roe_osc_pert(:, 4) * a_chief * 1e3)
plot(t_pert_1 / T, roe_mean_pert(:, 4) * a_chief * 1e3)
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(t_pert_1 / T, roe_osc_pert(:, 5) * a_chief * 1e3)
plot(t_pert_1 / T, roe_mean_pert(:, 5) * a_chief * 1e3)
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(t_pert_1 / T, roe_osc_pert(:, 6) * a_chief * 1e3)
plot(t_pert_1 / T, roe_mean_pert(:, 6) * a_chief * 1e3)
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

%% Question 4: Relative position in the RTN frame

rel_pos_unpert = zeros(N, 3);
rel_pos_pert = zeros(N, 3);

for i=1:N
    rel_state_ECI = y_unpert_1(i, 7:12) - y_unpert_1(i, 1:6);
    rel_state_RTN = ECI2RTN(y_unpert_1(i, 1:6), rel_state_ECI, mu);
    rel_pos_unpert(i, :) = rel_state_RTN(1:3)*1e3;

    rel_state_ECI = y_pert_1(i, 7:12) - y_pert_1(i, 1:6);
    rel_state_RTN = ECI2RTN(y_pert_1(i, 1:6), rel_state_ECI, mu);
    rel_pos_pert(i, :) = rel_state_RTN(1:3)*1e3;
end

figure
subplot(2,2,1)
hold on
plot(rel_pos_unpert(:, 2), rel_pos_unpert(:, 1), 'LineWidth', 1.5)
plot(rel_pos_pert(:, 2), rel_pos_pert(:, 1))
hold on
grid on
axis equal
legend('No J2', 'J2')
xlabel('T-axis [m]')
ylabel('R-axis [m]')

subplot(2,2,2)
hold on
plot(rel_pos_unpert(:, 3), rel_pos_unpert(:, 1), 'LineWidth', 1.5)
plot(rel_pos_pert(:, 3), rel_pos_pert(:, 1))
hold on
grid on
axis equal
xlabel('N-axis [m]')
ylabel('R-axis [m]')

subplot(2,2,3)
hold on
plot(rel_pos_unpert(:, 2), rel_pos_unpert(:, 3), 'LineWidth', 1.5)
plot(rel_pos_pert(:, 2), rel_pos_pert(:, 3))
hold on
grid on
axis equal
xlabel('T-axis [m]')
ylabel('N-axis [m]')

subplot(2,2,4)
hold on
plot3(rel_pos_unpert(:, 1), rel_pos_unpert(:, 2), rel_pos_unpert(:, 3), 'LineWidth', 1.5)
plot3(rel_pos_pert(:, 1), rel_pos_pert(:, 2), rel_pos_pert(:, 3))
hold off
grid on
axis equal
view(3)
xlabel('R-axis [m]')
ylabel('T-axis [m]')
zlabel('N-axis [km]')

%% Question 5: Planar plots of the ROE

figure
subplot(3,1,1)
hold on
plot(roe_osc_unpert(:,3) * a_chief * 1e3, roe_osc_unpert(:, 4) * a_chief * 1e3, 'LineWidth', 3)
plot(roe_mean_unpert(:, 3) * a_chief * 1e3, roe_mean_unpert(:, 4) * a_chief * 1e3)
% plot(0, 0,'*')
hold on
grid on
axis equal
legend('Osculatin', 'Mean')
xlabel('\deltae_x [m]')
ylabel('\deltae_y [m]')

subplot(3,1,2)
hold on
plot(roe_osc_unpert(:, 5) * a_chief * 1e3, roe_osc_unpert(:, 6) * a_chief * 1e3, 'LineWidth', 3)
plot(roe_mean_unpert(:, 5) * a_chief * 1e3, roe_mean_unpert(:, 6) * a_chief * 1e3)
hold off
grid on
axis equal
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')

subplot(3,1,3)
hold on
plot(roe_osc_unpert(:, 2) * a_chief * 1e3, roe_osc_unpert(:, 1) * a_chief * 1e3, 'LineWidth', 3)
plot(roe_mean_unpert(:, 2) * a_chief * 1e3, roe_mean_unpert(:, 1) * a_chief * 1e3)
hold off
grid on
axis equal
xlabel('\delta\lambda [m]')
ylabel('\deltaa [m]')

figure
subplot(3,1,1)
hold on
plot(roe_osc_pert(:,3) * a_chief * 1e3, roe_osc_pert(:, 4) * a_chief * 1e3)
plot(roe_mean_pert(:, 3) * a_chief * 1e3, roe_mean_pert(:, 4) * a_chief * 1e3)
hold on
grid on
axis equal
legend('Osculatin', 'Mean')
xlabel('\deltae_x [m]')
ylabel('\deltae_y [m]')

subplot(3,1,2)
hold on
plot(roe_osc_pert(:, 5) * a_chief * 1e3, roe_osc_pert(:, 6) * a_chief * 1e3)
plot(roe_mean_pert(:, 5) * a_chief * 1e3, roe_mean_pert(:, 6) * a_chief * 1e3)
hold off
grid on
axis equal
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')

subplot(3,1,3)
hold on
plot(roe_osc_pert(:, 2) * a_chief * 1e3, roe_osc_pert(:, 1) * a_chief * 1e3)
plot(roe_mean_pert(:, 2) * a_chief * 1e3, roe_mean_pert(:, 1) * a_chief * 1e3)
hold off
grid on
axis equal
xlabel('\delta\lambda [m]')
ylabel('\deltaa [m]')

%% Question 7: Initial Conditions after Maneuver

init_mean_oe_deputy_2 = [initial_mean_oe_chief(1), initial_mean_oe_deputy(2), initial_mean_oe_chief(3), initial_mean_oe_deputy(5), initial_mean_oe_deputy(4), initial_mean_oe_deputy(6)];
init_osc_oe_deputy_2 = mean2osc(init_mean_oe_deputy_2);

% Computing the size of the maneuvers
(init_osc_oe_deputy_2(1) - a_deputy*1e3) * n_deputy / 2
(init_osc_oe_deputy_2(3) - inc_deputy) * n_deputy * a_deputy * 1e3

delta_a = (init_osc_oe_deputy_2(1) * 1e-3 - a_chief) / a_chief;
delta_lambda = (init_osc_oe_deputy_2(6) + init_osc_oe_deputy_2(5)) - (M_chief + omega_chief) + (init_osc_oe_deputy_2(4) - RAAN_chief) * cos(inc_chief);
delta_ex = init_osc_oe_deputy_2(2) * cos(init_osc_oe_deputy_2(5)) - e_chief * cos(omega_chief);
delta_ey = init_osc_oe_deputy_2(2) * sin(init_osc_oe_deputy_2(5)) - e_chief * sin(omega_chief);
delta_ix = init_osc_oe_deputy_2(3) - inc_chief;
delta_iy = (init_osc_oe_deputy_2(4) - RAAN_chief) * sin(inc_chief); 

delta_a_mean = (init_mean_oe_deputy_2(1) - init_mean_oe_chief(1)) / init_mean_oe_chief(1);
delta_lambda_mean = (init_mean_oe_deputy_2(6) + init_mean_oe_deputy_2(5)) - (init_mean_oe_chief(6) + init_mean_oe_chief(5)) + (init_mean_oe_deputy_2(4) - init_mean_oe_chief(4)) * cos(init_mean_oe_chief(3));
delta_ex_mean = init_mean_oe_deputy_2(2) * cos(init_mean_oe_deputy_2(5)) - init_mean_oe_chief(2) * cos(init_mean_oe_chief(5));
delta_ey_mean = init_mean_oe_deputy_2(2) * sin(init_mean_oe_deputy_2(5)) - init_mean_oe_chief(2) * sin(init_mean_oe_chief(5));
delta_ix_mean = init_mean_oe_deputy_2(3) - init_mean_oe_chief(3);
delta_iy_mean = (init_mean_oe_deputy_2(4) - init_mean_oe_chief(4)) * sin(init_mean_oe_chief(3));

a_deputy = a_chief + delta_a * a_chief;
inc_deputy = inc_chief + delta_ix;
RAAN_deputy = RAAN_chief + delta_iy / sin(inc_chief);
e_deputy = sqrt((delta_ex + e_chief * cos(omega_chief))^2 + (delta_ey + e_chief * sin(omega_chief))^2);
omega_deputy = atan2(delta_ey + e_chief * sin(omega_chief), delta_ex + e_chief * cos(omega_chief));
M_deputy = delta_lambda + M_chief + omega_chief - (RAAN_deputy - RAAN_chief) * cos(inc_chief) - omega_deputy;
E_deputy = M2E(M_deputy, e_deputy, 1e-12);
nu_deputy = E2f(E_deputy, e_deputy);
init_oe_deputy_2 = [a_deputy * 1e3, e_deputy, inc_deputy, RAAN_deputy, omega_deputy, M_deputy];

[pos_deputy_2, vel_deputy_2] = OE2ECI(a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, nu_deputy, mu);

init_cond_2 = [pos_chief', vel_chief', pos_deputy_2', vel_deputy_2'];

% Numerical integration without J2 perturbations
[t_unpert_2, y_unpert_2] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, init_cond_2, options);

roe_osc_unpert_2 = zeros(N, 6);

roe_mean_unpert_2 = zeros(N, 6);

for i=1:N
    [a_c, e_c, i_c, omega_c, RAAN_c, M_c] = ECI2OE_M(y_unpert_2(i,1:6), mu);
    oe_c = [a_c, e_c, i_c, omega_c, RAAN_c, M_c];
    [a_d, e_d, i_d, omega_d, RAAN_d, M_d] = ECI2OE_M(y_unpert_2(i,7:12), mu);
    oe_d = [a_d, e_d, i_d, omega_d, RAAN_d, M_d];
    ROE = OE2ROE(oe_c, oe_d);
    roe_osc_unpert_2(i, :) = ROE;

    oe_c_mean = osc2mean([oe_c(1)*1e3, oe_c(2:3), oe_c(5), oe_c(4), oe_c(6)], 0);
    oe_d_mean = osc2mean([oe_d(1)*1e3, oe_d(2:3), oe_d(5), oe_d(4), oe_d(6)], 0);
    oe_c_mean = [oe_c_mean(1)*1e-3, oe_c_mean(2:3)', oe_c_mean(5), oe_c_mean(4), oe_c_mean(6)];
    oe_d_mean = [oe_d_mean(1)*1e-3, oe_d_mean(2:3)', oe_d_mean(5), oe_d_mean(4), oe_d_mean(6)];
    ROE_mean = OE2ROE(oe_c_mean, oe_d_mean);
    roe_mean_unpert_2(i, :) = ROE_mean;
end

figure
subplot(3,2,1)
hold on
plot(t_unpert_1 / T, roe_osc_unpert_2(:, 1) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert_2(:, 1) * a_chief * 1e3, 'r--')
hold off
grid on
legend('Osculating', 'Mean')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(t_unpert_1 / T, roe_osc_unpert_2(:, 2) * a_chief * 1e3, 'b-', 'LineWidth', 2)
plot(t_unpert_1 / T, roe_mean_unpert_2(:, 2) * a_chief * 1e3, 'r--', 'LineWidth', 2)
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(t_unpert_1 / T, roe_osc_unpert_2(:, 3) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert_2(:, 3) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_unpert_1 / T, roe_osc_unpert_2(:, 4) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert_2(:, 4) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(t_unpert_1 / T, roe_osc_unpert_2(:, 5) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert_2(:, 5) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(t_unpert_1 / T, roe_osc_unpert_2(:, 6) * a_chief * 1e3, 'b-')
plot(t_unpert_1 / T, roe_mean_unpert_2(:, 6) * a_chief * 1e3, 'r--')
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(roe_osc_unpert_2(:,3) * a_chief * 1e3, roe_osc_unpert_2(:, 4) * a_chief * 1e3, 'LineWidth', 2)
plot(roe_mean_unpert_2(:, 3) * a_chief * 1e3, roe_mean_unpert_2(:, 4) * a_chief * 1e3,'--', 'LineWidth', 2)
hold on
grid on
axis equal
legend('Osculatin', 'Mean')
xlabel('\deltae_x [m]')
ylabel('\deltae_y [m]')

subplot(3,1,2)
hold on
plot(roe_osc_unpert_2(:, 5) * a_chief * 1e3, roe_osc_unpert_2(:, 6) * a_chief * 1e3, 'LineWidth', 2)
plot(roe_mean_unpert_2(:, 5) * a_chief * 1e3, roe_mean_unpert_2(:, 6) * a_chief * 1e3,'--', 'LineWidth', 2)
hold off
grid on
axis equal
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')

subplot(3,1,3)
hold on
plot(roe_osc_unpert_2(:, 2) * a_chief * 1e3, roe_osc_unpert_2(:, 1) * a_chief * 1e3, 'LineWidth', 2)
plot(roe_mean_unpert_2(:, 2) * a_chief * 1e3, roe_mean_unpert_2(:, 1) * a_chief * 1e3,'--', 'LineWidth', 2)
hold off
grid on
axis equal
xlabel('\delta\lambda [m]')
ylabel('\deltaa [m]')

% Numerical integration with J2 perturbations
[t_pert_2, y_pert_2] = ode89(@(t, state) state_der_posvel_J2_2sats(t, state, mu, J2, R_E), tspan, init_cond_2, options);

roe_osc_pert_2 = zeros(N, 6);

roe_mean_pert_2 = zeros(N, 6);

for i=1:N
    [a_c, e_c, i_c, omega_c, RAAN_c, M_c] = ECI2OE_M(y_pert_2(i,1:6), mu);
    oe_c = [a_c, e_c, i_c, omega_c, RAAN_c, M_c];
    [a_d, e_d, i_d, omega_d, RAAN_d, M_d] = ECI2OE_M(y_pert_2(i,7:12), mu);
    oe_d = [a_d, e_d, i_d, omega_d, RAAN_d, M_d];
    ROE = OE2ROE(oe_c, oe_d);
    roe_osc_pert_2(i, :) = ROE;

    oe_c_mean = osc2mean([oe_c(1)*1e3, oe_c(2:3), oe_c(5), oe_c(4), oe_c(6)]);
    oe_d_mean = osc2mean([oe_d(1)*1e3, oe_d(2:3), oe_d(5), oe_d(4), oe_d(6)]);
    oe_c_mean = [oe_c_mean(1)*1e-3, oe_c_mean(2:3)', oe_c_mean(5), oe_c_mean(4), oe_c_mean(6)];
    oe_d_mean = [oe_d_mean(1)*1e-3, oe_d_mean(2:3)', oe_d_mean(5), oe_d_mean(4), oe_d_mean(6)];
    ROE_mean = OE2ROE(oe_c_mean, oe_d_mean);
    roe_mean_pert_2(i, :) = ROE_mean;
end

figure
subplot(3,2,1)
hold on
plot(t_pert_2 / T, roe_osc_pert_2(:, 1) * a_chief * 1e3)
plot(t_pert_2 / T, roe_mean_pert_2(:, 1) * a_chief * 1e3)
hold off
grid on
legend('Osculating', 'Mean')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(t_pert_2 / T, roe_osc_pert_2(:, 2) * a_chief * 1e3)
plot(t_pert_2 / T, roe_mean_pert_2(:, 2) * a_chief * 1e3)
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(t_pert_2 / T, roe_osc_pert_2(:, 3) * a_chief * 1e3)
plot(t_pert_2 / T, roe_mean_pert_2(:, 3) * a_chief * 1e3)
hold off
grid on
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_pert_2 / T, roe_osc_pert_2(:, 4) * a_chief * 1e3)
plot(t_pert_2 / T, roe_mean_pert_2(:, 4) * a_chief * 1e3)
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(t_pert_2 / T, roe_osc_pert_2(:, 5) * a_chief * 1e3)
plot(t_pert_2 / T, roe_mean_pert_2(:, 5) * a_chief * 1e3)
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(t_pert_2 / T, roe_osc_pert_2(:, 6) * a_chief * 1e3)
plot(t_pert_2 / T, roe_mean_pert_2(:, 6) * a_chief * 1e3)
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(roe_osc_pert_2(:,3) * a_chief * 1e3, roe_osc_pert_2(:, 4) * a_chief * 1e3)
plot(roe_mean_pert_2(:, 3) * a_chief * 1e3, roe_mean_pert_2(:, 4) * a_chief * 1e3, 'LineWidth', 2)
hold on
grid on
axis equal
legend('Osculatin', 'Mean')
xlabel('\deltae_x [m]')
ylabel('\deltae_y [m]')

subplot(3,1,2)
hold on
plot(roe_osc_pert_2(:, 5) * a_chief * 1e3, roe_osc_pert_2(:, 6) * a_chief * 1e3)
plot(roe_mean_pert_2(:, 5) * a_chief * 1e3, roe_mean_pert_2(:, 6) * a_chief * 1e3,'.', 'LineWidth', 2)
hold off
grid on
axis equal
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')

subplot(3,1,3)
hold on
plot(roe_osc_pert_2(:, 2) * a_chief * 1e3, roe_osc_pert_2(:, 1) * a_chief * 1e3)
plot(roe_mean_pert_2(:, 2) * a_chief * 1e3, roe_mean_pert_2(:, 1) * a_chief * 1e3, 'LineWidth', 2)
hold off
grid on
axis equal
xlabel('\delta\lambda [m]')
ylabel('\deltaa [m]')

%% Question 8: Using the STM in both scenarios

% Computing the initial mean ROE before applying the maneuver
init_mean_roe = OE2ROE(initial_mean_oe_chief, initial_mean_oe_deputy);

roe_bef_man = zeros(N,6);

for i=1:N
    e = sqrt(oe_c_mean_pert(i, 3)^2 + oe_c_mean_pert(i, 4)^2);
    omega = atan2(oe_c_mean_pert(i,4), oe_c_mean_pert(i, 3));
    STM = STM_QNSROE(tspan(i), oe_c_mean_pert(i, 1), e, oe_c_mean_pert(i, 5), omega, J2, R_E, mu);
    roe_bef_man(i,:) = STM * init_mean_roe';
end

figure
subplot(3,1,1)
hold on
plot(roe_osc_pert(:,3) * a_chief * 1e3, roe_osc_pert(:, 4) * a_chief * 1e3)
plot(roe_mean_pert(:, 3) * a_chief * 1e3, roe_mean_pert(:, 4) * a_chief * 1e3, 'LineWidth', 2)
plot(roe_bef_man(:, 3) * a_chief * 1e3, roe_bef_man(:, 4) * a_chief * 1e3, '--', 'LineWidth', 2)
hold on
grid on
axis equal
legend('Osculatin', 'Mean', 'STM')
xlabel('\deltae_x [m]')
ylabel('\deltae_y [m]')

subplot(3,1,2)
hold on
plot(roe_osc_pert(:, 5) * a_chief * 1e3, roe_osc_pert(:, 6) * a_chief * 1e3)
plot(roe_mean_pert(:, 5) * a_chief * 1e3, roe_mean_pert(:, 6) * a_chief * 1e3, 'LineWidth', 2)
plot(roe_bef_man(:, 5) * a_chief * 1e3, roe_bef_man(:, 6) * a_chief * 1e3,'--','LineWidth', 2)
hold off
grid on
axis equal
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')

subplot(3,1,3)
hold on
plot(roe_osc_pert(:, 2) * a_chief * 1e3, roe_osc_pert(:, 1) * a_chief * 1e3)
plot(roe_mean_pert(:, 2) * a_chief * 1e3, roe_mean_pert(:, 1) * a_chief * 1e3, 'LineWidth', 2)
plot(roe_bef_man(:, 2) * a_chief * 1e3, roe_bef_man(:, 1) * a_chief * 1e3,'--','LineWidth', 2)
hold off
grid on
axis equal
xlabel('\delta\lambda [m]')
ylabel('\deltaa [m]')

figure
subplot(3,2,1)
plot(t_unpert_2 / T, abs(roe_mean_pert(:, 1) - roe_bef_man(:, 1)) * a_chief * 1e3)
ylabel('\deltaa [m]')

subplot(3,2,3)
plot(t_unpert_2 / T, abs(roe_mean_pert(:, 2) - roe_bef_man(:, 2)) * a_chief * 1e3)
ylabel('\delta\lambda [m]')

subplot(3,2,5)
plot(t_unpert_2 / T, abs(roe_mean_pert(:, 3) - roe_bef_man(:, 3)) * a_chief * 1e3)
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
plot(t_unpert_2 / T, abs(roe_mean_pert(:, 4) - roe_bef_man(:, 4)) * a_chief * 1e3)
ylabel('\deltae_y [m]')

subplot(3,2,4)
plot(t_unpert_2 / T, abs(roe_mean_pert(:, 5) - roe_bef_man(:, 5)) * a_chief * 1e3)
ylabel('\deltai_x [m]')

subplot(3,2,6)
plot(t_unpert_2 / T, abs(roe_mean_pert(:, 6) - roe_bef_man(:, 6)) * a_chief * 1e3)
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')


% Computing the initial mean ROE after applying the maneuver
init_mean_oe_deputy = osc2mean(init_oe_deputy_2);
initial_mean_oe_deputy = [init_mean_oe_deputy(1:3)', init_mean_oe_deputy(5), init_mean_oe_deputy(4), init_mean_oe_deputy(6)];
init_mean_roe_2 = OE2ROE(initial_mean_oe_chief, initial_mean_oe_deputy);
% init_mean_roe = OE2ROE(initial_mean_oe_chief, [init_mean_oe_deputy_2(1:3), init_mean_oe_deputy_2(5), init_mean_oe_deputy_2(4), init_mean_oe_deputy_2(6)]);

roe_aft_man = zeros(N,6);

for i=1:N
    e = sqrt(oe_c_mean_pert(i, 3)^2 + oe_c_mean_pert(i, 4)^2);
    omega = atan2(oe_c_mean_pert(i,4), oe_c_mean_pert(i, 3));
    STM = STM_QNSROE(tspan(i), oe_c_mean_pert(i, 1), e, oe_c_mean_pert(i, 5), omega, J2, R_E, mu);
    roe_aft_man(i,:) = STM * init_mean_roe_2';
end

figure
subplot(3,1,1)
hold on
plot(roe_osc_pert_2(:,3) * a_chief * 1e3, roe_osc_pert_2(:, 4) * a_chief * 1e3)
plot(roe_mean_pert_2(:, 3) * a_chief * 1e3, roe_mean_pert_2(:, 4) * a_chief * 1e3, 'LineWidth', 2)
plot(roe_aft_man(:, 3) * a_chief * 1e3, roe_aft_man(:, 4) * a_chief * 1e3,'--', 'LineWidth', 2)
hold on
grid on
axis equal
legend('Osculatin', 'Mean', 'STM')
xlabel('\deltae_x [m]')
ylabel('\deltae_y [m]')

subplot(3,1,2)
hold on
plot(roe_osc_pert_2(:, 5) * a_chief * 1e3, roe_osc_pert_2(:, 6) * a_chief * 1e3)
plot(roe_mean_pert_2(:, 5) * a_chief * 1e3, roe_mean_pert_2(:, 6) * a_chief * 1e3,'o', 'LineWidth', 4)
plot(roe_aft_man(:, 5) * a_chief * 1e3, roe_aft_man(:, 6) * a_chief * 1e3, 'o', 'LineWidth', 2)
hold off
grid on
axis equal
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')

subplot(3,1,3)
hold on
plot(roe_osc_pert_2(:, 2) * a_chief * 1e3, roe_osc_pert_2(:, 1) * a_chief * 1e3)
plot(roe_mean_pert_2(:, 2) * a_chief * 1e3, roe_mean_pert_2(:, 1) * a_chief * 1e3, 'LineWidth', 3)
plot(roe_aft_man(:, 2) * a_chief * 1e3, roe_aft_man(:, 1) * a_chief * 1e3,'o', 'LineWidth', 3)
hold off
grid on
axis equal
xlabel('\delta\lambda [m]')
ylabel('\deltaa [m]')

figure
subplot(3,2,1)
plot(t_unpert_2 / T, abs(roe_mean_pert_2(:, 1) - roe_aft_man(:, 1)) * a_chief * 1e3)
ylabel('\deltaa [m]')

subplot(3,2,3)
plot(t_unpert_2 / T, abs(roe_mean_pert_2(:, 2) - roe_aft_man(:, 2)) * a_chief * 1e3)
ylabel('\delta\lambda [m]')

subplot(3,2,5)
plot(t_unpert_2 / T, abs(roe_mean_pert_2(:, 3) - roe_aft_man(:, 3)) * a_chief * 1e3)
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
plot(t_unpert_2 / T, abs(roe_mean_pert_2(:, 4) - roe_aft_man(:, 4)) * a_chief * 1e3)
ylabel('\deltae_y [m]')

subplot(3,2,4)
plot(t_unpert_2 / T, abs(roe_mean_pert_2(:, 5) - roe_aft_man(:, 5)) * a_chief * 1e3)
ylabel('\deltai_x [m]')

subplot(3,2,6)
plot(t_unpert_2 / T, abs(roe_mean_pert_2(:, 6) - roe_aft_man(:, 6)) * a_chief * 1e3)
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

%% Functions

function statedot = state_der_posvel_J2_2sats(t, state, mu, J2, R_E)
    statedot = zeros(size(state));
    statedot(1:3) = state(4:6);
    statedot(7:9) = state(10:12);
    r1 = norm(state(1:3));
    d = - mu .* J2 .* R_E.^2 ./ 2 .* (6 .* state(3) ./ r1.^5 .* [0 0 1] + (3 ./ r1.^4 - 15 .* state(3).^2 ./ r1.^6) .* state(1:3)' ./ r1);
    statedot(4:6) = - mu .* state(1:3) ./ r1.^3 + d';
    r2 = norm(state(7:9));
    d = - mu .* J2 .* R_E.^2 ./ 2 .* (6 .* state(9) ./ r2.^5 .* [0 0 1] + (3 ./ r2.^4 - 15 .* state(9).^2 ./ r2.^6) .* state(7:9)' ./ r2);
    statedot(10:12) = - mu .* state(7:9) ./ r2.^3 + d';
end

function statedot = FODE_2sats(t, state, mu)
    statedot = zeros(size(state));
    statedot(1:3) = state(4:6);
    statedot(7:9) = state(10:12);

    r0 = state(1:3);
    r1 = state(7:9);

    statedot(4:6) = - mu * r0 / norm(r0)^3;
    statedot(10:12) = - mu * r1 / norm(r1)^3;
end

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapTo2Pi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end

function QNSOE = OE2QNSOE(oe)
    QNSOE = zeros(1, 6);
    QNSOE(1) = oe(1);
    QNSOE(2) = wrapTo2Pi(oe(4) + oe(6));
    QNSOE(3) = oe(2) * cos(oe(4));
    QNSOE(4) = oe(2) * sin(oe(4));
    QNSOE(5) = oe(3);
    QNSOE(6) = oe(5);
end

function STM = STM_QNSROE(t, a, e, i, omega, J2, R_E, mu)
    STM = eye(6);
    eta = sqrt(1 - e^2);
    n = sqrt(mu / a^3);
    k = 3 * J2 * R_E^2 * sqrt(mu) / (4 * a^(7/2) * eta^4);
    E = 1 + eta;
    F = 4 + 3 * eta;
    G = 1 / eta^2;
    P = 3 * cos(i)^2 - 1;
    Q = 5 * cos(i)^2 - 1;
    S = sin(2*i);
    T = sin(i)^2;

    e_xi = e * cos(omega);
    e_yi = e * sin(omega);
    omega_dot = k * Q;
    omega_f = omega + omega_dot * t;
    e_xf = e * cos(omega_f);
    e_yf = e * sin(omega_f);

    STM(2,1) = -(3/2 * n + 7/2 * k * E * P) * t;
    STM(2,3) = k * e_xi * F * G * P * t;
    STM(2,4) = k * e_yi * F * G * P * t;
    STM(2,5) = -k * F * S * t;

    STM(3,1) = 7/2 * k * e_yf * Q * t;
    STM(3,3) = cos(omega_dot * t) - 4 * k * e_xi * e_yf * G * Q * t;
    STM(3,4) = - sin(omega_dot * t) - 4 *k * e_yi * e_yf * G * Q * t;
    STM(3,5) = 5 * k * e_yf * S * t;

    STM(4,1) = -7/2 * k * e_xf * Q * t;
    STM(4,3) = sin(omega_dot * t) + 4 * k * e_xi * e_xf * G * Q * t;
    STM(4,4) = cos(omega_dot * t) + 4 * k * e_yi * e_xf * G * Q * t;
    STM(4,5) = -5 *k * e_xf * S * t;

    STM(6,1) = 7/2 * k * S * t;
    STM(6,3) = -4 * k * e_xi * G * S * t;
    STM(6,4) = -4 * k * e_yi * G * S * t;
    STM(6,5) = 2 * k * T * t;
end