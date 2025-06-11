%% Problem 1

% This is trying with the sensor giving the ROE straight away

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 1); % options for numerical integration

% Initial conditions = the ones at the start of the perigee pass mode
a_chief = 36943; % km
e_chief = 0.8111;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
T_chief = 2 * pi / n_chief;
M_chief = n_chief * (14 * 3600 + 49 * 60); % perigee pass mode starts at t=14h49
nu_chief_PPM = mean2true(M_chief, e_chief, 1e-12);
oe_chief_PPM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

delta_a = 0;
delta_lambda = 100e-3 / a_chief;
delta_ex = 750e-3 / a_chief;
delta_ey = 150e-3 / a_chief;
delta_ix = 700e-3 / a_chief;
delta_iy = 140e-3 / a_chief;

a_deputy = a_chief + delta_a * a_chief;
inc_deputy = inc_chief + delta_ix;
RAAN_deputy = RAAN_chief + delta_iy / sin(inc_chief);
e_deputy = sqrt((delta_ex + e_chief * cos(omega_chief))^2 + (delta_ey + e_chief * sin(omega_chief))^2);
omega_deputy = atan2(delta_ey + e_chief * sin(omega_chief), delta_ex + e_chief * cos(omega_chief));
M_deputy = delta_lambda + M_chief + omega_chief - (RAAN_deputy - RAAN_chief) * cos(inc_chief) - omega_deputy;
nu_deputy_PPM = mean2true(M_deputy, e_deputy, 1e-12);
oe_deputy_PPM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

ROE_PPM = OE2ROE(oe_chief_PPM, oe_deputy_PPM);

% Defining the noise matrices
% P_0 = zeros(6);
% P_0(1,1) = (5e-8 / (3 * a_chief))^2;
% P_0(2,2) = (5e-7 / (3 * a_chief))^2;
% P_0(3,3) = (1e-9 / (3 * a_chief))^2;
% P_0(4,4) = (1e-9 / (3 * a_chief))^2;
% P_0(5,5) = (7e-10 / (3 * a_chief))^2;
% P_0(6,6) = (5e-10 / (3 * a_chief))^2;
% P_0_sqrt = sqrtm(P_0);
% 
% Q = zeros(6);
% Q(1,1) = (3e-9 / (3 * a_chief))^2;
% Q(2,2) = (3e-8 / (3 * a_chief))^2;
% Q(3,3) = (6e-10 / (3 * a_chief))^2;
% Q(4,4) = (7e-10 / (3 * a_chief))^2;
% Q(5,5) = (3e-10 / (3 * a_chief))^2;
% Q(6,6) = (1e-10 / (3 * a_chief))^2;
% 
% R = eye(6) * (1e-9 / (3 * a_chief))^2;
% R(1,1) = (5e-9 / (3 * a_chief))^2;
% R(2,2) = (5e-8 / (3 * a_chief))^2;
% R(5,5) = (5e-10 / (3 * a_chief))^2;
% R(6,6) = (1e-10 / (3 * a_chief))^2;

P_0 = zeros(6);
P_0(1,1) = (5e-9 / (3 * a_chief))^2;
P_0(2,2) = (5e-8 / (3 * a_chief))^2;
P_0(3,3) = (1e-9 / (3 * a_chief))^2;
P_0(4,4) = (1e-7 / (3 * a_chief))^2;
P_0(5,5) = (3e-9 / (3 * a_chief))^2;
P_0(6,6) = (5e-8 / (3 * a_chief))^2;
P_0 = 1e2 * P_0;
P_0_sqrt = sqrtm(P_0);

% Q = zeros(6);
% Q(1,1) = (3e-9 / (3 * a_chief))^2;
% Q(2,2) = (3e-8 / (3 * a_chief))^2;
% Q(3,3) = (6e-10 / (3 * a_chief))^2;
% Q(4,4) = (5e-9 / (3 * a_chief))^2;
% Q(5,5) = (5e-9 / (3 * a_chief))^2;
% Q(6,6) = (5e-10 / (3 * a_chief))^2;
% Q = 1e-2 * P_0;
Q = zeros(6);
Q(1,1) = (3e-9 / (3 * a_chief))^2;
Q(2,2) = (3e-8 / (3 * a_chief))^2;
Q(3,3) = (6e-10 / (3 * a_chief))^2;
Q(4,4) = (5e-7 / (3 * a_chief))^2;
Q(5,5) = (5e-9 / (3 * a_chief))^2;
Q(6,6) = (5e-8 / (3 * a_chief))^2;

R = eye(6) * (5e-10 / (3 * a_chief))^2;
R(1,1) = (1e-8 / (3 * a_chief))^2;
R(2,2) = (1e-7 / (3 * a_chief))^2;
R(4,4) = (2e-7 / (3 * a_chief))^2;
R(5,5) = (4e-9 / (3 * a_chief))^2;
R(6,6) = (5e-8 / (3 * a_chief))^2;
R = R * 1e1;

% Initial conditions (for the filter and ground truth)
ic = sqrtm(P_0)*randn(6,1) + ROE_PPM';
[pos_chief, vel_chief] = OE2ECI(oe_chief_PPM(1), oe_chief_PPM(2), oe_chief_PPM(3), oe_chief_PPM(4), oe_chief_PPM(5), nu_chief_PPM, mu);
[pos_deputy, vel_deputy] = OE2ECI(oe_deputy_PPM(1), oe_deputy_PPM(2), oe_deputy_PPM(3), oe_deputy_PPM(4), oe_deputy_PPM(5), nu_deputy_PPM, mu);
ic_ECI = [pos_chief', vel_chief', pos_deputy', vel_deputy'];

rel_pos_gt = pos_deputy - pos_chief;
rel_vel_gt = vel_deputy - vel_chief;
rel_posvel_gt_RTN = ECI2RTN([pos_chief', vel_chief'], [rel_pos_gt', rel_vel_gt'], mu);

oe_deputy_PPM_est = ROE2OE(oe_chief_PPM, ic);
nu_deputy_PPM_est = mean2true(oe_deputy_PPM_est(6), oe_deputy_PPM_est(2));
[pos_deputy_est, vel_deputy_est] = OE2ECI(oe_deputy_PPM_est(1), oe_deputy_PPM_est(2), oe_deputy_PPM_est(3), oe_deputy_PPM_est(4), oe_deputy_PPM_est(5), nu_deputy_PPM_est, mu);

rel_pos_est = pos_deputy_est - pos_chief;
rel_vel_est = vel_deputy_est - vel_chief;
rel_posvel_est_RTN = ECI2RTN([rel_pos_est', rel_vel_est'], [pos_chief', vel_chief'], mu);

% Simulation of the ground truth
N = 1000; % 200000
tspan = linspace(0, 1 * T_chief, N); % 20
[t_ECI, y_ECI] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, ic_ECI, options);

% Filter
history_ROE_filter = zeros(N, 6);
history_ROE_filter(1, :) = ic;
history_ROE_gt = zeros(N, 6);
history_ROE_gt(1, :) = ROE_PPM;
history_covariance = zeros(6, 6, N);
history_covariance(:, :, 1) = P_0;
history_sqrt_covariance = zeros(N, 6); % for just the diagonal entries
history_sqrt_covariance(1, 1) = P_0_sqrt(1,1);
history_sqrt_covariance(1, 2) = P_0_sqrt(2,2);
history_sqrt_covariance(1, 3) = P_0_sqrt(3,3);
history_sqrt_covariance(1, 4) = P_0_sqrt(4,4);
history_sqrt_covariance(1, 5) = P_0_sqrt(5,5);
history_sqrt_covariance(1, 6) = P_0_sqrt(6,6);
history_prefit = zeros(N-1, 6);
history_postfit = zeros(N-1, 6);
history_rel_posvel_gt = zeros(N, 6);
history_rel_posvel_gt(1, :) = rel_posvel_gt_RTN;
history_rel_posvel_est = zeros(N, 6);
history_rel_posvel_est(1, :) = rel_posvel_est_RTN;

for j=2:N
    j
    % propagate the state and covariance (no control inputs for now)
    prev = history_ROE_filter(j-1, :);
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j-1, 1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j-1, 7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = OE2ROE(oe_chief, oe_deputy);
    rel_pos_gt = y_ECI(j-1, 7:9) - y_ECI(j-1, 1:3);
%     norm(rel_pos_gt)
    rel_vel_gt = y_ECI(j-1, 10:12) - y_ECI(j-1, 4:6);
    rel_posvel_gt_RTN = ECI2RTN(y_ECI(j-1, 1:6), [rel_pos_gt, rel_vel_gt], mu);
    history_ROE_gt(j, :) = ROE;
    history_rel_posvel_gt(j-1, :) = rel_posvel_gt_RTN;
    oe_deputy_est = ROE2OE(oe_chief, prev);
    nu_deputy_est = mean2true(oe_deputy_est(6), oe_deputy_est(2), 1e-12);
    [pos_deputy_est, vel_deputy_est] = OE2ECI(oe_deputy_est(1), oe_deputy_est(2), oe_deputy_est(3), oe_deputy_est(4), oe_deputy_est(5), nu_deputy_est, mu);
    
    rel_pos_est = pos_deputy_est' - y_ECI(j-1, 1:3);
    rel_vel_est = vel_deputy_est' - y_ECI(j-1, 4:6);
    rel_posvel_est_RTN = ECI2RTN(y_ECI(j-1, 1:6), [rel_pos_est, rel_vel_est], mu);
    history_rel_posvel_est(j-1, :) = rel_posvel_est_RTN;
    span = linspace(tspan(j-1), tspan(j), 2);

    % integrating to propagate the state
    [t, y] = ode89(@(t, state) ROE_noJ2(t, state, mu, a_chief, oe_deputy_est(1)), span, prev, options);
    x_pred = y(end, :) + (sqrtm(Q) * randn(6, 1))';

    % using the STM to propagate the covariance
    A = ROE_STM(tspan(j)-tspan(j-1), n_chief);
    P = A * history_covariance(:, :, j-1) * A' + Q;

    % refit with the measurements (H is just identity for this one so
    % excluded from the formula)
    ROE_meas = sqrtm(R)*randn(6,1) + ROE';
    K = P * (P + R)^(-1);
    new_state_est = x_pred' + K * (ROE_meas - x_pred');
    history_ROE_filter(j, :) = new_state_est;
    new_P = (eye(6) - K) * P * (eye(6) - K)' + K * R * K';
    history_covariance(:, :, j) = new_P;
    new_P_sqrt = sqrtm(new_P);
    history_sqrt_covariance(j, 1) = new_P_sqrt(1,1);
    history_sqrt_covariance(j, 2) = new_P_sqrt(2,2);
    history_sqrt_covariance(j, 3) = new_P_sqrt(3,3);
    history_sqrt_covariance(j, 4) = new_P_sqrt(4,4);
    history_sqrt_covariance(j, 5) = new_P_sqrt(5,5);
    history_sqrt_covariance(j, 6) = new_P_sqrt(6,6);

    % Prefit residual
    pre_res = ROE_meas - x_pred';
    history_prefit(j-1, :) = pre_res;

    % Postfit residual
    post_res = ROE_meas - new_state_est;
    history_postfit(j-1, :) = post_res;
end

%% Plots

figure
subplot(3,2,1)
hold on
plot(tspan / T_chief, history_ROE_filter(:, 1) * a_chief * 1e3)
plot(tspan / T_chief, history_ROE_gt(:, 1) * a_chief * 1e3)
hold off
grid on
legend('Filter estimation', 'Ground truth')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(tspan / T_chief, history_ROE_filter(:, 2) * a_chief * 1e3)
plot(tspan / T_chief, history_ROE_gt(:, 2) * a_chief * 1e3)
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(tspan / T_chief, history_ROE_filter(:, 3) * a_chief * 1e3)
plot(tspan / T_chief, history_ROE_gt(:, 3) * a_chief * 1e3)
hold off
grid on
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(tspan / T_chief, history_ROE_filter(:, 4) * a_chief * 1e3)
plot(tspan / T_chief, history_ROE_gt(:, 4) * a_chief * 1e3)
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(tspan / T_chief, history_ROE_filter(:, 5) * a_chief * 1e3)
plot(tspan / T_chief, history_ROE_gt(:, 5) * a_chief * 1e3)
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(tspan / T_chief, history_ROE_filter(:, 6) * a_chief * 1e3)
plot(tspan / T_chief, history_ROE_gt(:, 6) * a_chief * 1e3)
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

% 
% figure
% subplot(3,2,1)
% plot(tspan / T_chief, (history_ROE_gt(:, 1) - history_ROE_filter(:, 1)) * a_chief * 1e3)
% grid on
% ylabel('\deltaa [m]')
% 
% subplot(3,2,3)
% plot(tspan / T_chief, (history_ROE_gt(:, 2) - history_ROE_filter(:, 2)) * a_chief * 1e3)
% grid on
% ylabel('\delta\lambda [m]')
% 
% subplot(3,2,5)
% plot(tspan / T_chief, (history_ROE_gt(:, 3) - history_ROE_filter(:, 3)) * a_chief * 1e3)
% grid on
% ylabel('\deltae_x [m]')
% xlabel('Orbital Periods')
% 
% subplot(3,2,2)
% plot(tspan / T_chief, (history_ROE_gt(:, 4) - history_ROE_filter(:, 4)) * a_chief * 1e3)
% grid on
% ylabel('\deltae_y [m]')
% 
% subplot(3,2,4)
% plot(tspan / T_chief, (history_ROE_gt(:, 5) - history_ROE_filter(:, 5)) * a_chief * 1e3)
% grid on
% ylabel('\deltai_x [m]')
% 
% subplot(3,2,6)
% plot(tspan / T_chief, (history_ROE_gt(:, 6) - history_ROE_filter(:, 6)) * a_chief * 1e3)
% grid on
% ylabel('\deltai_y [m]')
% xlabel('Orbital Periods')


figure
subplot(3,2,1)
hold on
plot(tspan / T_chief, (history_ROE_gt(:, 1) - history_ROE_filter(:, 1)) * a_chief * 1e3)
plot(tspan / T_chief, 3 * history_sqrt_covariance(:, 1) * a_chief * 1e3, 'r-')
plot(tspan / T_chief, -3 * history_sqrt_covariance(:, 1) * a_chief * 1e3, 'r-')
hold off
legend('Estimation error', 'Covariance bound (3\sigma)')
grid on
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(tspan / T_chief, (history_ROE_gt(:, 2) - history_ROE_filter(:, 2)) * a_chief * 1e3)
plot(tspan / T_chief, 3 * history_sqrt_covariance(:, 2) * a_chief * 1e3, 'r-')
plot(tspan / T_chief, -3 * history_sqrt_covariance(:, 2) * a_chief * 1e3, 'r-')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(tspan / T_chief, (history_ROE_gt(:, 3) - history_ROE_filter(:, 3)) * a_chief * 1e3)
plot(tspan / T_chief, 3 * history_sqrt_covariance(:, 3) * a_chief * 1e3, 'r-')
plot(tspan / T_chief, -3 * history_sqrt_covariance(:, 3) * a_chief * 1e3, 'r-')
hold off
grid on
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(tspan / T_chief, (history_ROE_gt(:, 4) - history_ROE_filter(:, 4)) * a_chief * 1e3)
plot(tspan / T_chief, 3 * history_sqrt_covariance(:, 4) * a_chief * 1e3, 'r-')
plot(tspan / T_chief, -3 * history_sqrt_covariance(:, 4) * a_chief * 1e3, 'r-')
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(tspan / T_chief, (history_ROE_gt(:, 5) - history_ROE_filter(:, 5)) * a_chief * 1e3)
plot(tspan / T_chief, 3 * history_sqrt_covariance(:, 5) * a_chief * 1e3, 'r-')
plot(tspan / T_chief, -3 * history_sqrt_covariance(:, 5) * a_chief * 1e3, 'r-')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(tspan / T_chief, (history_ROE_gt(:, 6) - history_ROE_filter(:, 6)) * a_chief * 1e3)
plot(tspan / T_chief, 3 * history_sqrt_covariance(:, 6) * a_chief * 1e3, 'r-')
plot(tspan / T_chief, -3 * history_sqrt_covariance(:, 6) * a_chief * 1e3, 'r-')
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')


figure
subplot(3,2,1)
hold on
plot(tspan(2:end) / T_chief, history_prefit(:, 1) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, history_postfit(:, 1) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, ones(N-1, 1) * sqrt(R(1,1)) * 3 * a_chief * 1e3, 'k-')
plot(tspan(2:end) / T_chief, -ones(N-1, 1) * sqrt(R(1,1)) * 3 * a_chief * 1e3, 'k-')
hold off
grid on 
legend('Prefit', 'Postfit', 'Measurement noise (3\sigma)')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(tspan(2:end) / T_chief, history_prefit(:, 2) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, history_postfit(:, 2) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, ones(N-1, 1) * sqrt(R(2,2)) * 3 * a_chief * 1e3, 'k-')
plot(tspan(2:end) / T_chief, -ones(N-1, 1) * sqrt(R(2,2)) * 3 * a_chief * 1e3, 'k-')
hold off
grid on 
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(tspan(2:end) / T_chief, history_prefit(:, 3) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, history_postfit(:, 3) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, ones(N-1, 1) * sqrt(R(3,3)) * 3 * a_chief * 1e3, 'k-')
plot(tspan(2:end) / T_chief, -ones(N-1, 1) * sqrt(R(3,3)) * 3 * a_chief * 1e3, 'k-')
hold off
grid on 
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(tspan(2:end) / T_chief, history_prefit(:, 4) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, history_postfit(:, 4) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, ones(N-1, 1) * sqrt(R(4,4)) * 3 * a_chief * 1e3, 'k-')
plot(tspan(2:end) / T_chief, -ones(N-1, 1) * sqrt(R(4,4)) * 3 * a_chief * 1e3, 'k-')
hold off
grid on 
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(tspan(2:end) / T_chief, history_prefit(:, 5) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, history_postfit(:, 5) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, ones(N-1, 1) * sqrt(R(5,5)) * 3 * a_chief * 1e3, 'k-')
plot(tspan(2:end) / T_chief, -ones(N-1, 1) * sqrt(R(5,5)) * 3 * a_chief * 1e3, 'k-')
hold off
grid on 
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(tspan(2:end) / T_chief, history_prefit(:, 6) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, history_postfit(:, 6) * a_chief * 1e3)
plot(tspan(2:end) / T_chief, ones(N-1, 1) * sqrt(R(6,6)) * 3 * a_chief * 1e3, 'k-')
plot(tspan(2:end) / T_chief, -ones(N-1, 1) * sqrt(R(6,6)) * 3 * a_chief * 1e3, 'k-')
hold off
grid on 
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')


figure
subplot(3,2,1)
plot(tspan / T_chief, (history_rel_posvel_gt(:, 1) - history_rel_posvel_est(:, 1)) * 1e3)
grid on
ylabel('Error along R [m]')

subplot(3,2,3)
plot(tspan / T_chief, (history_rel_posvel_gt(:, 2) - history_rel_posvel_est(:, 2)) * 1e3)
grid on
ylabel('Error along T [m]')

subplot(3,2,5)
plot(tspan / T_chief, (history_rel_posvel_gt(:, 3) - history_rel_posvel_est(:, 3)) * 1e3)
grid on
ylabel('Error along N [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
plot(tspan / T_chief, (history_rel_posvel_gt(:, 4) - history_rel_posvel_est(:, 4)) * 1e3)
grid on
ylabel('Error along R [m/s]')

subplot(3,2,4)
plot(tspan / T_chief, (history_rel_posvel_gt(:, 5) - history_rel_posvel_est(:, 5)) * 1e3)
grid on
ylabel('Error along T [m/s]')

subplot(3,2,6)
plot(tspan / T_chief, (history_rel_posvel_gt(:, 6) - history_rel_posvel_est(:, 6)) * 1e3)
grid on
ylabel('Error along N [m/s]')
xlabel('Orbital Periods')


%% Functions

function statedot = ROE_noJ2(t, state, mu, a_c, a_d)
    statedot = zeros(size(state));
    statedot(2) = sqrt(mu) * (1/sqrt(a_d^3) - 1/sqrt(a_c^3));
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
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end

function oe_deputy = ROE2OE(oe_chief, ROE)
    oe_deputy = zeros(6,1);
    oe_deputy(1) = oe_chief(1) + ROE(1) * oe_chief(1);
    oe_deputy(3) = wrapTo2Pi(oe_chief(3) + ROE(5));
    oe_deputy(2) = sqrt((ROE(3)+oe_chief(2)*cos(oe_chief(4)))^2 + (ROE(4)+oe_chief(2)*sin(oe_chief(4)))^2);
    oe_deputy(5) = wrapTo2Pi(ROE(6) / sin(oe_chief(3)) + oe_chief(5));
    oe_deputy(4) = wrapTo2Pi(atan2(ROE(4)+oe_chief(2)*sin(oe_chief(4)), ROE(3)+oe_chief(2)*cos(oe_chief(4))));
    oe_deputy(6) = wrapTo2Pi(ROE(2) + oe_chief(6) - oe_deputy(4) + oe_chief(4) - (oe_deputy(5) - oe_chief(5))*cos(oe_chief(3)));
end

function A = ROE_STM(t, n)
    A = eye(6);
    A(2,1) = -1.5 * n * t;
end