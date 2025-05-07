%% For station keeping

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
R_E = 6378; % km
J2 = 0.108263e-2;

%% Problem 1

% Defining the ROE for Inertial Attitude Mode (IAM)

% Initial conditions in mean elements for the chief

a_chief = 36943; % km
e_chief = 0.8111;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
T = 2 * pi / n_chief;
M_chief = n_chief * (6 * 3600 + 49 * 60);
nu_chief_IAM = mean2true(M_chief, e_chief);
oe_mean_chief_IAM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];
temp = mean2osc([oe_mean_chief_IAM(1)*1e3 oe_mean_chief_IAM(2:3) oe_mean_chief_IAM(5) oe_mean_chief_IAM(4) oe_mean_chief_IAM(6)]);
oe_osc_chief_IAM = [temp(1)*1e-3 temp(2:3)' temp(5) temp(4) temp(6)];
nu_osc_chief_IAM = mean2true(temp(6), temp(2));

% Initial conditions in mean elements for the deputy

delta_a = 0;
delta_lambda = 80e-3 / a_chief;
delta_ex = -20e-3 / a_chief;
delta_ey = 89.4e-3 / a_chief;
delta_ix = 0 / a_chief;
delta_iy = 79e-3 / a_chief;

a_deputy = a_chief + delta_a * a_chief;
inc_deputy = inc_chief + delta_ix;
RAAN_deputy = RAAN_chief + delta_iy / sin(inc_chief);
e_deputy = sqrt((delta_ex + e_chief * cos(omega_chief))^2 + (delta_ey + e_chief * sin(omega_chief))^2);
omega_deputy = atan2(delta_ey + e_chief * sin(omega_chief), delta_ex + e_chief * cos(omega_chief));
M_deputy = delta_lambda + M_chief + omega_chief - (RAAN_deputy - RAAN_chief) * cos(inc_chief) - omega_deputy;
nu_deputy_IAM = mean2true(M_deputy, e_deputy);
oe_mean_deputy_IAM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];
temp = mean2osc([oe_mean_deputy_IAM(1)*1e3 oe_mean_deputy_IAM(2:3) oe_mean_deputy_IAM(5) oe_mean_deputy_IAM(4) oe_mean_deputy_IAM(6)]);
oe_osc_deputy_IAM = [temp(1)*1e-3 temp(2:3)' temp(5) temp(4) temp(6)];
nu_osc_deputy_IAM = mean2true(temp(6), temp(2));

init_ROE_IAM = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]';
init_mod_ROE_IAM = OE2ROE_modified(oe_mean_chief_IAM, oe_mean_deputy_IAM)';

% Checking initial separation and initial conditions for FODE
[pos_chief, vel_chief] = OE2ECI(oe_osc_chief_IAM(1), oe_osc_chief_IAM(2), oe_osc_chief_IAM(3), oe_osc_chief_IAM(4), oe_osc_chief_IAM(5), nu_osc_chief_IAM, mu);
[pos_deputy, vel_deputy] = OE2ECI(oe_osc_deputy_IAM(1), oe_osc_deputy_IAM(2), oe_osc_deputy_IAM(3), oe_osc_deputy_IAM(4), oe_osc_deputy_IAM(5), nu_osc_deputy_IAM, mu);
rho_ECI = norm(pos_deputy - pos_chief);
norm(pos_chief);

% Numerical integration without J2 perturbations
init_FODE = [pos_chief, vel_chief, pos_deputy, vel_deputy];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100);
N = 100000;
tspan = linspace(0, 100 * T, N);
[t, y] = ode89(@(t, state) FODE_2sats_J2(t, state, mu, J2, R_E), tspan, init_FODE, options);

mod_ROE_gt = zeros(N,6);
mod_ROE_stm = zeros(N,6);
rel_motion_RTN_gt = zeros(N,6);
rel_motion_RTN_stm = zeros(N,6);

for j=1:N
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y(j,1:6), mu);
    oe_osc_chief = [a, e, i, omega, RAAN, M];
    temp = osc2mean([oe_osc_chief(1)*1e3 oe_osc_chief(2:3) oe_osc_chief(5) oe_osc_chief(4) oe_osc_chief(6)]);
    oe_mean_chief = [temp(1)*1e-3 temp(2:3)' temp(5) temp(4) temp(6)];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y(j,7:12), mu);
    oe_osc_deputy = [a, e, i, omega, RAAN, M];
    temp = osc2mean([oe_osc_deputy(1)*1e3 oe_osc_deputy(2:3) oe_osc_deputy(5) oe_osc_deputy(4) oe_osc_deputy(6)]);
    oe_mean_deputy = [temp(1)*1e-3 temp(2:3)' temp(5) temp(4) temp(6)];
    mod_ROE = OE2ROE_modified(oe_mean_chief, oe_mean_deputy);
    mod_ROE_gt(j,:) = mod_ROE;

%     rel_state_ECI = y(j,7:12) - y(j,1:6);
%     rel_state_RTN = ECI2RTN(y(j,1:6), rel_state_ECI, mu);
    nu_mean_chief_cur = mean2true(oe_mean_chief(6), oe_mean_chief(2));
    nu_mean_deputy_cur = mean2true(oe_mean_deputy(6), oe_mean_deputy(2));
    [pos_c, vel_c] = OE2ECI(oe_mean_chief(1), oe_mean_chief(2), oe_mean_chief(3), oe_mean_chief(4), oe_mean_chief(5), nu_mean_chief_cur, mu);
    [pos_d, vel_d] = OE2ECI(oe_mean_deputy(1), oe_mean_deputy(2), oe_mean_deputy(3), oe_mean_deputy(4), oe_mean_deputy(5), nu_mean_deputy_cur, mu);
    rel_state_ECI = [pos_d', vel_d'] - [pos_c', vel_c'];
    rel_state_RTN = ECI2RTN([pos_c', vel_c'], rel_state_ECI, mu);
    rel_motion_RTN_gt(j, :) = rel_state_RTN;

    STM = STM_ROE_J2(tspan(j), oe_mean_chief(1), oe_mean_chief(2), oe_mean_chief(3), oe_mean_chief_IAM(4), J2, R_E, mu);
    ROE = STM * init_ROE_IAM;
    oe_mean_deputy = ROE2OE(oe_mean_chief, ROE);
    mod_ROE = OE2ROE_modified(oe_mean_chief, oe_mean_deputy);
    nu_mean_chief_cur = mean2true(oe_mean_chief(6), oe_mean_chief(2));
    nu_mean_deputy_cur = mean2true(oe_mean_deputy(6), oe_mean_deputy(2));
    [pos_c, vel_c] = OE2ECI(oe_mean_chief(1), oe_mean_chief(2), oe_mean_chief(3), oe_mean_chief(4), oe_mean_chief(5), nu_mean_chief_cur, mu);
    [pos_d, vel_d] = OE2ECI(oe_mean_deputy(1), oe_mean_deputy(2), oe_mean_deputy(3), oe_mean_deputy(4), oe_mean_deputy(5), nu_mean_deputy_cur, mu);
    rel_state_ECI = [pos_d', vel_d'] - [pos_c', vel_c'];
    rel_state_RTN = ECI2RTN([pos_c', vel_c'], rel_state_ECI, mu);
    rel_motion_RTN_stm(j, :) = rel_state_RTN;
    mod_ROE_stm(j, :) = mod_ROE;
end

figure
subplot(3,2,1)
hold on
plot(t / T, mod_ROE_gt(:, 1) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 1) * a_chief * 1e3, '--')
hold off
grid on
legend('Ground Truth', 'STM')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(t / T, mod_ROE_gt(:, 2) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 2) * a_chief * 1e3, '--')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(t / T, mod_ROE_gt(:, 3) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 3) * a_chief * 1e3, '--')
hold off
grid on
ylabel("\deltae_x' [m]")
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t / T, mod_ROE_gt(:, 4) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 4) * a_chief * 1e3, '--')
hold off
grid on
ylabel("\deltae_y' [m]")

subplot(3,2,4)
hold on
plot(t / T, mod_ROE_gt(:, 5) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 5) * a_chief * 1e3, '--')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(t / T, mod_ROE_gt(:, 6) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 6) * a_chief * 1e3, '--')
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

figure
subplot(3,1,1)
hold on
plot(mod_ROE_gt(:, 2) * a_chief * 1e3, mod_ROE_gt(:, 1) * a_chief * 1e3)
plot(mod_ROE_stm(:, 2) * a_chief * 1e3, mod_ROE_stm(:, 1) * a_chief * 1e3, '--')
hold off
axis equal
grid on
legend('Ground Truth', 'STM')
xlabel('\delta\lambda [m]')
ylabel('\deltaa [m]')

subplot(3,1,2)
hold on
plot(mod_ROE_gt(:,3) * a_chief * 1e3, mod_ROE_gt(:,4) * a_chief * 1e3)
plot(mod_ROE_stm(:,3) * a_chief * 1e3, mod_ROE_stm(:,4) * a_chief * 1e3, '--')
hold off
axis equal
grid on
xlabel("\deltae_x' [m]")
ylabel("\deltae_y' [m]")

subplot(3,1,3)
hold on
plot(mod_ROE_gt(:,5) * a_chief * 1e3, mod_ROE_gt(:,6) * a_chief * 1e3)
plot(mod_ROE_stm(:,5) * a_chief * 1e3, mod_ROE_stm(:,6) * a_chief * 1e3, '--')
hold off
axis equal
grid on
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')

figure
subplot(2,2,1)
plot(rel_motion_RTN_gt(:,2), rel_motion_RTN_gt(:,1))
axis equal
grid on
xlabel('T-axis [km]')
ylabel('R-axis [km]')

subplot(2,2,2)
plot(rel_motion_RTN_gt(:,3), rel_motion_RTN_gt(:,1))
axis equal
grid on
xlabel('N-axis [km]')
ylabel('R-axis [km]')

subplot(2,2,3)
plot(rel_motion_RTN_gt(:,2), rel_motion_RTN_gt(:,3))
axis equal
grid on
xlabel('T-axis [km]')
ylabel('N-axis [km]')

subplot(2,2,4)
plot3(rel_motion_RTN_gt(:,1), rel_motion_RTN_gt(:,2), rel_motion_RTN_gt(:,3))
view(3)
axis equal
grid on
xlabel('R-axis [km]')
ylabel('T-axis [km]')
zlabel('N-axis [km]')

figure
subplot(2,2,1)
plot(rel_motion_RTN_gt(:,5), rel_motion_RTN_gt(:,4))
axis equal
grid on
xlabel('T-axis [km/s]')
ylabel('R-axis [km/s]')

subplot(2,2,2)
plot(rel_motion_RTN_gt(:,6), rel_motion_RTN_gt(:,4))
axis equal
grid on
xlabel('N-axis [km/s]')
ylabel('R-axis [km/s]')

subplot(2,2,3)
plot(rel_motion_RTN_gt(:,5), rel_motion_RTN_gt(:,6))
axis equal
grid on
xlabel('T-axis [km/s]')
ylabel('N-axis [km/s]')

subplot(2,2,4)
plot3(rel_motion_RTN_gt(:,4), rel_motion_RTN_gt(:,5), rel_motion_RTN_gt(:,6))
view(3)
axis equal
grid on
xlabel('R-axis [km/s]')
ylabel('T-axis [km/s]')
zlabel('N-axis [km/S]')

figure
subplot(3,2,1)
plot(t / T, abs(rel_motion_RTN_gt(:,1) - rel_motion_RTN_stm(:,1)))
ylabel('Eror in R-axis [km]')

subplot(3,2,3)
plot(t / T, abs(rel_motion_RTN_gt(:,2) - rel_motion_RTN_stm(:,2)))
ylabel('Error in T-axis [km]')

subplot(3,2,5)
plot(t / T, abs(rel_motion_RTN_gt(:,3) - rel_motion_RTN_stm(:,3)))
ylabel('Error in N-axis [km]')
xlabel("Orbital Periods")

subplot(3,2,2)
plot(t / T, abs(rel_motion_RTN_gt(:,4) - rel_motion_RTN_stm(:,4)))
ylabel('Error in R-axis [km/s]')

subplot(3,2,4)
plot(t / T, abs(rel_motion_RTN_gt(:,5) - rel_motion_RTN_stm(:,5)))
ylabel('Error in T-axis [km/s]')

subplot(3,2,6)
plot(t / T, abs(rel_motion_RTN_gt(:,6) - rel_motion_RTN_stm(:,6)))
ylabel('Error in N-axis [km/s]')
xlabel('Orbital Periods')

%% Problem 2

%% Functions

function statedot = FODE_2sats_J2(t, state, mu, J2, R_E)
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

function STM = STM_ROE_J2(t, a, e, i, omega, J2, R_E, mu)
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

function ROE = OE2ROE_modified(oe_chief, oe_deputy)
%     eta = sqrt(1 - oe_chief(2)^2);
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
%     ROE(2) = wrapToPi((oe_deputy(6) + eta*oe_deputy(4)) - (oe_chief(6) + eta*oe_chief(4)) + eta*(oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) - oe_chief(2);
    ROE(4) = wrapToPi(oe_deputy(4) - oe_chief(4) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = wrapToPi((oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3)));
end

function oe_deputy = ROE2OE(oe_chief, ROE)
    oe_deputy = zeros(1,6);
    oe_deputy(1) = oe_chief(1) + oe_chief(1) * ROE(1);
    oe_deputy(2) = sqrt((ROE(3)+oe_chief(2)*cos(oe_chief(4)))^2 + (ROE(4) + oe_chief(2) * sin(oe_chief(4)))^2);
    oe_deputy(4) = atan2(ROE(4)+oe_chief(2)*sin(oe_chief(4)), ROE(3)+oe_chief(2)*cos(oe_chief(4)));
    oe_deputy(3) = oe_chief(3) + ROE(5);
    oe_deputy(5) = oe_chief(5) + ROE(6) / sin(oe_chief(3));
    oe_deputy(6) = ROE(2) - oe_deputy(4) + oe_chief(6) + oe_chief(4) - (oe_deputy(5) - oe_chief(5))*cos(oe_chief(3));
end