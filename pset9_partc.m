%% Problem 2

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
M_chief = n_chief * (4 * 3600 + 49 * 60); % perigee pass mode starts at t=14h49
nu_chief_PPM = mean2true(M_chief, e_chief, 1e-15);
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
nu_deputy_PPM = mean2true(M_deputy, e_deputy, 1e-15);
oe_deputy_PPM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

ROE_PPM = OE2ROE(oe_chief_PPM, oe_deputy_PPM);

% % Computing the final modified ROE for the reconfiguration (Inertial
% % Attitude Mode)

a_chief = 36943; % km
e_chief = 0.8111;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
T = 2 * pi / n_chief;
M_chief = n_chief * (6 * 3600 + 49 * 60); % reconfiguration ends at t=6h49
% nu_chief_IAM = mean2true(M_chief, e_chief);
oe_chief_IAM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

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
% nu_deputy_IAM = mean2true(M_deputy, e_deputy);
oe_deputy_IAM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

modified_ROE_IAM = eccentric_ROE(oe_chief_IAM, oe_deputy_IAM);
ROE_IAM = OE2ROE(oe_chief_IAM, oe_deputy_IAM);

% Defining the noise matrices
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
ic = sqrtm(P_0) * randn(6,1) + ROE_PPM';
[pos_chief, vel_chief] = OE2ECI(oe_chief_PPM(1), oe_chief_PPM(2), oe_chief_PPM(3), oe_chief_PPM(4), oe_chief_PPM(5), nu_chief_PPM, mu);
[pos_deputy, vel_deputy] = OE2ECI(oe_deputy_PPM(1), oe_deputy_PPM(2), oe_deputy_PPM(3), oe_deputy_PPM(4), oe_deputy_PPM(5), nu_deputy_PPM, mu);
ic_ECI = [pos_chief', vel_chief', pos_deputy', vel_deputy'];

rel_pos_gt = pos_deputy - pos_chief;
rel_vel_gt = vel_deputy - vel_chief;

oe_deputy_PPM_est = ROE2OE(oe_chief_PPM, ic);
nu_deputy_PPM_est = mean2true(oe_deputy_PPM_est(6), oe_deputy_PPM_est(2));
[pos_deputy_est, vel_deputy_est] = OE2ECI(oe_deputy_PPM_est(1), oe_deputy_PPM_est(2), oe_deputy_PPM_est(3), oe_deputy_PPM_est(4), oe_deputy_PPM_est(5), nu_deputy_PPM_est, mu);

rel_pos_est = pos_deputy_est - pos_chief;
rel_vel_est = vel_deputy_est - vel_chief;
% rel_posvel_est_RTN = ECI2RTN([rel_pos_est', rel_vel_est'], [pos_chief', vel_chief'], mu);

% Simulation of the ground truth
N = 10000;
tspan = linspace(0, 0.1 * T_chief, N);
% [t_ECI, y_ECI] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, ic_ECI, options);

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
history_rel_posvel_gt(1, :) = [rel_pos_gt', rel_vel_gt'];
history_rel_posvel_est = zeros(N, 6);
history_rel_posvel_est(1, :) = [rel_pos_est', rel_vel_est'];
control_tracking_error_history = zeros(N, 6);
control_tracking_error_history(1, :) = ROE_PPM - ROE_IAM;
deltav_history = zeros(N, 2);
history_posvel_chief = zeros(N, 6);
history_posvel_chief(1, :) = [pos_chief', vel_chief'];
history_posvel_deputy = zeros(N, 6);
history_posvel_deputy(1, :) = [pos_deputy', vel_deputy'];

prev = ic;

for j=2:N
    j
    % propagate the state and covariance
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j-1, 1:6), mu);
    [a, e, i, omega, RAAN, M] = ECI2OE_M(history_posvel_chief(j-1,:), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    nu_chief = mean2true(M, e, 1e-14);

    oe_deputy_est = ROE2OE(oe_chief, prev);
%     oe_deputy_est(2)
%     oe_deputy_est(6)
    nu_deputy_est = mean2true(oe_deputy_est(6), oe_deputy_est(2), 1e-14);

    % Compute the control law
    modified_ROE_cur = eccentric_ROE(oe_chief, oe_deputy_est);
    delta_alpha = modified_ROE_cur - modified_ROE_IAM;
    delta_alpha_nonmod = prev - ROE_IAM';
    control_tracking_error_history(j,:) = delta_alpha_nonmod;

    % If not controlling delta_lambda, applied = goal
    modified_ROE_applied = modified_ROE_IAM;

    % Controlling delta_lambda
    if abs(delta_alpha(2)) > 1e-6
        % From Lippe
        delta_a_applied = sign(delta_alpha(2)) * min([abs(delta_alpha(2)) / (500e-3 / oe_chief(1)), 100e-3 / oe_chief(1)]);
        modified_ROE_applied(1) = delta_a_applied;
    end

    delta_alpha_applied = modified_ROE_cur - modified_ROE_applied;

    B_mod = modified_control_matrix(oe_chief(1), sqrt(mu / oe_chief(1)^3), oe_chief(2), nu_chief, oe_chief(4));

    % Computing the optimal location of the out-of-plane maneuver
    nu_oop = wrapToPi(atan2(delta_alpha(6), delta_alpha(5)) - oe_chief(4));

    % Computing the optimal location for the in-plane maneuver
    delta_ex_tild = delta_alpha_applied(3);
    delta_ey_tild = oe_chief(2) * delta_alpha_applied(4);
    if delta_ex_tild == 0
        nu_ip = acos((sqrt(1-oe_chief(2)^2) - 1) / oe_chief(2));
    elseif delta_ey_tild == 0
        nu_ip = 0;
    else
        del_ex = oe_deputy_est(2) * cos(oe_deputy_est(4)) - oe_chief(2) * cos(oe_chief(4)) - (oe_deputy_IAM(2) * cos(oe_deputy_IAM(4)) - oe_chief_IAM(2) * cos(oe_chief_IAM(4)));
        del_ey = oe_deputy_est(2) * sin(oe_deputy_est(4)) - oe_chief(2) * sin(oe_chief(4)) - (oe_deputy_IAM(2) * sin(oe_deputy_IAM(4)) - oe_chief_IAM(2) * sin(oe_chief_IAM(4)));
        nu_ip = cast(wrapToPi(inplane_planning(oe_chief(2), oe_chief(4), oe_chief(3), del_ex, del_ey, delta_alpha_applied(6)) - oe_chief(4)), 'double');
    end

    G = control_gain(3000, 4, nu_deputy_est, nu_ip, nu_oop);
    u = - pinv(B_mod) * (zeros(5) * [modified_ROE_cur(1); modified_ROE_cur(3:end)'] + G * [delta_alpha_applied(1); delta_alpha_applied(3:end)']);
    u = u .* (tspan(j) - tspan(j-1)); % in km/s
    % should be 2D because no radial maneuvers 
    deltav_history(j+1,:) = u; % in km/s

    span = linspace(tspan(j-1), tspan(j), 2);

    % integrating to propagate the estimated state
    [t, y] = ode89(@(t, state) ROE_noJ2_control(t, state, mu, oe_chief, u / (tspan(j) - tspan(j-1))), span, prev, options);
%     [t, y] = ode89(@(t, state) ROE_noJ2_control(t, state, mu, oe_chief, [0, 0]), span, prev, options);
%     B = control_matrix(oe_chief(1), sqrt(mu / oe_chief(1)^3), oe_chief(2), nu_chief, oe_chief(4), oe_chief(3));
%     B = control_matrix(oe_deputy_est(1), sqrt(mu / oe_deputy_est(1)^3), oe_deputy_est(2), nu_deputy_est, oe_deputy_est(4), oe_deputy_est(3));
    x_pred = y(end, :) + (sqrtm(Q) * randn(6, 1))'; % + (B * u)';

    % using the STM to propagate the covariance
    A = ROE_STM(tspan(j)-tspan(j-1), n_chief);
    P = A * history_covariance(:, :, j-1) * A' + Q; % nothing changes when we add control inputs?

    % Propagate the ground truth with the newly computed control input
    % Rotate the delta v from RTN to ECI (using deputy's RTN)
    % RTN unit vectors
    R_vec = history_posvel_deputy(j-1, 1:3) / norm(history_posvel_deputy(j-1, 1:3));
    h = cross(history_posvel_deputy(j-1, 1:3), history_posvel_deputy(j-1, 4:6));
    N_vec = h / norm(h);
    T = cross(N_vec, R_vec);

    % Rotation matrix from ECI to RTN
    rotation = [R_vec; T; N_vec];
    delta_v = rotation' * [0; u / (tspan(j) - tspan(j-1))];

%     tspan = [reconfiguration_tspan(j), reconfiguration_tspan(j+1)];
    init_cond = [history_posvel_chief(j-1, :), history_posvel_deputy(j-1,:)];
    [t, y] = ode89(@(t, state) FODE_2sats(t, state, mu, delta_v), span, init_cond, options);
%     [t, y] = ode89(@(t, state) FODE_2sats(t, state, mu, [0; 0; 0]), span, init_cond, options);
    history_posvel_chief(j, :) = y(end, 1:6);
    history_posvel_deputy(j, :) = y(end, 7:12);

    % refit with the measurements (H is just identity for this one so
    % excluded from the formula)
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j, 1:6), mu);
    [a, e, i, omega, RAAN, M] = ECI2OE_M(history_posvel_chief(j,:), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j, 7:12), mu);
    [a, e, i, omega, RAAN, M] = ECI2OE_M(history_posvel_deputy(j,:), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = OE2ROE(oe_chief, oe_deputy);
    ROE_meas = sqrtm(R)*randn(6,1) + ROE';
    K = P * (P + R)^(-1);
    new_state_est = x_pred' + K * (ROE_meas - x_pred');
    prev = new_state_est;
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

%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j, 1:6), mu);
%     oe_chief = [a, e, i, omega, RAAN, M];
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j, 7:12), mu);
%     oe_deputy = [a, e, i, omega, RAAN, M];
%     ROE = OE2ROE(oe_chief, oe_deputy);
    rel_pos_gt = history_posvel_deputy(j, 1:3) - history_posvel_chief(j, 1:3);
    rel_vel_gt = history_posvel_deputy(j, 4:6) - history_posvel_chief(j, 4:6);
%     rel_posvel_gt_RTN = ECI2RTN([rel_pos_gt, rel_vel_gt], y_ECI(j-1, 1:6), mu);
    history_ROE_gt(j, :) = ROE;
%     history_rel_posvel_gt(j-1, :) = rel_posvel_gt_RTN;
    history_rel_posvel_gt(j, :) = [rel_pos_gt, rel_vel_gt];

    oe_deputy_est = ROE2OE(oe_chief, new_state_est);
    nu_deputy_est = mean2true(oe_deputy_est(6), oe_deputy_est(2), 1e-14);
    [pos_deputy_est, vel_deputy_est] = OE2ECI(oe_deputy_est(1), oe_deputy_est(2), oe_deputy_est(3), oe_deputy_est(4), oe_deputy_est(5), nu_deputy_est, mu);
    
    rel_pos_est = pos_deputy_est' - history_posvel_chief(j, 1:3);
    rel_vel_est = vel_deputy_est' - history_posvel_chief(j, 4:6);
%     rel_posvel_est_RTN = ECI2RTN([rel_pos_est, rel_vel_est], y_ECI(j-1, 1:6), mu);
%     history_rel_posvel_est(j-1, :) = rel_posvel_est_RTN;
    history_rel_posvel_est(j, :) = [rel_pos_est, rel_vel_est];
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


% Don't forget that this is not in the RTN frame but in ECI for now
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

figure
subplot(3,2,1)
hold on
plot(tspan(1) / T_chief, ROE_PPM(1) * a_chief * 1e3, 'ro')
plot(tspan / T_chief, control_tracking_error_history(:,1) * a_chief * 1e3, 'b-')
plot(tspan(end) / T_chief, ROE_IAM(1) * a_chief * 1e3, 'rx')
hold off
grid on
legend('Start', 'Maneuvering', 'Target')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(tspan(1) / T_chief, ROE_PPM(2) * a_chief * 1e3, 'ro')
plot(tspan / T_chief, control_tracking_error_history(:,2) * a_chief * 1e3, 'b-')
plot(tspan(end) / T_chief, ROE_IAM(2) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(tspan(1) / T_chief, ROE_PPM(3) * a_chief * 1e3, 'ro')
plot(tspan / T_chief, control_tracking_error_history(:,3) * a_chief * 1e3, 'b-')
plot(tspan(end) / T_chief, ROE_IAM(3) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltae_x [m]')

subplot(3,2,2)
hold on
plot(tspan(1) / T_chief, ROE_PPM(4) * a_chief * 1e3, 'ro')
plot(tspan / T_chief, control_tracking_error_history(:,4) * a_chief * 1e3, 'b-')
plot(tspan(end) / T_chief, ROE_IAM(4) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(tspan(1) / T_chief, ROE_PPM(5) * a_chief * 1e3, 'ro')
plot(tspan / T_chief, control_tracking_error_history(:,5) * a_chief * 1e3, 'b-')
plot(tspan(end) / T_chief, ROE_IAM(5) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(tspan(1) / T_chief, ROE_PPM(6) * a_chief * 1e3, 'ro')
plot(tspan / T_chief, control_tracking_error_history(:,6) * a_chief * 1e3, 'b-')
plot(tspan(end) / T_chief, ROE_IAM(6) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltai_y [m]')

%% Functions

function statedot = ROE_noJ2(t, state, mu, a_c, a_d)
    statedot = zeros(size(state));
    statedot(2) = sqrt(mu) * (1/sqrt(a_d^3) - 1/sqrt(a_c^3));
end

function statedot = ROE_noJ2_control(t, state, mu, oe_chief, u)
    oe_deputy_estimated = ROE2OE(oe_chief, state);
    a_d = oe_deputy_estimated(1);
    n_d = sqrt(mu / a_d^3);
    e_d = oe_deputy_estimated(2);
    i_d = oe_deputy_estimated(3);
    omega_d = oe_deputy_estimated(4);
    RAAN_d = oe_deputy_estimated(5);
    M_d = oe_deputy_estimated(6);
    nu_d = mean2true(M_d, e_d, 1e-14);
    r_d = a_d * (1 - e_d^2) / (1 + e_d * cos(nu_d));

    a_c = oe_chief(1);
    n_c = sqrt(mu / oe_chief(1)^3);
    i_c = oe_chief(3);

    % u should be 2-dimensional

    statedot = zeros(size(state));
    
%     statedot(1) = 2 * a_d * sqrt(1 - e_d^2) / (n_d * r_d) * u(1);
%     statedot(2) = (n_d - (1 - e_d^2) / (n_d * a_d * e_d) * (1 + r_d / (a_d * (1 - e_d^2))) * sin(nu_d) * u(1) + (sqrt(1 - e_d^2) / (n_d * a_d * e_d) * (2 + e_d * cos(nu_d)) / (1 + e_d * cos(nu_d)) * sin(nu_d) * u(1) - r_d * cot(i_d) * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2)) * u(2)) - n_c + cos(i_c) * r_d * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2) * sin(i_d)) * u(2)) * a_c;
%     statedot(3) = (sqrt(1 - e_d^2) / (n_d * a_d^2 * e_d) * (a_d^2 * (1 - e_d^2) / r_d - r_d) * cos(omega_d) * u(1) - e_d * (sqrt(1 - e_d^2) / (n_d * a_d * e_d) * (2 + e_d * cos(nu_d)) / (1 + e_d * cos(nu_d)) * sin(nu_d) * u(1) - r_d * cot(i_d) * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2)) * u(2)) * sin(omega_d)) * a_c;
%     statedot(4) = (sqrt(1 - e_d^2) / (n_d * a_d^2 * e_d) * (a_d^2 * (1 - e_d^2) / r_d - r_d) * sin(omega_d) * u(1) + e_d * (sqrt(1 - e_d^2) / (n_d * a_d * e_d) * (2 + e_d * cos(nu_d)) / (1 + e_d * cos(nu_d)) * sin(nu_d) * u(1) - r_d * cot(i_d) * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2)) * u(2)) * cos(omega_d)) * a_c;
%     statedot(5) = a_c * r_d * cos(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2)) * u(2);
%     statedot(6) = a_c * sin(i_c) * r_d * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2) * sin(i_d)) * u(2);
    statedot(1) = 2 * a_d * sqrt(1 - e_d^2) / (n_d * r_d) * u(1) / a_c;
    statedot(2) = (n_d - (1 - e_d^2) / (n_d * a_d * e_d) * (1 + r_d / (a_d * (1 - e_d^2))) * sin(nu_d) * u(1) + (sqrt(1 - e_d^2) / (n_d * a_d * e_d) * (2 + e_d * cos(nu_d)) / (1 + e_d * cos(nu_d)) * sin(nu_d) * u(1) - r_d * cot(i_d) * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2)) * u(2)) - n_c + cos(i_c) * r_d * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2) * sin(i_d)) * u(2));
    statedot(3) = (sqrt(1 - e_d^2) / (n_d * a_d^2 * e_d) * (a_d^2 * (1 - e_d^2) / r_d - r_d) * cos(omega_d) * u(1) - e_d * (sqrt(1 - e_d^2) / (n_d * a_d * e_d) * (2 + e_d * cos(nu_d)) / (1 + e_d * cos(nu_d)) * sin(nu_d) * u(1) - r_d * cot(i_d) * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2)) * u(2)) * sin(omega_d));
    statedot(4) = (sqrt(1 - e_d^2) / (n_d * a_d^2 * e_d) * (a_d^2 * (1 - e_d^2) / r_d - r_d) * sin(omega_d) * u(1) + e_d * (sqrt(1 - e_d^2) / (n_d * a_d * e_d) * (2 + e_d * cos(nu_d)) / (1 + e_d * cos(nu_d)) * sin(nu_d) * u(1) - r_d * cot(i_d) * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2)) * u(2)) * cos(omega_d));
    statedot(5) = r_d * cos(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2)) * u(2);
    statedot(6) = sin(i_c) * r_d * sin(omega_d + nu_d) / (n_d * a_d^2 * sqrt(1 - e_d^2) * sin(i_d)) * u(2);
end

% function statedot = FODE_2sats(t, state, mu)
%     statedot = zeros(size(state));
%     statedot(1:3) = state(4:6);
%     statedot(7:9) = state(10:12);
% 
%     r0 = state(1:3);
%     r1 = state(7:9);
% 
%     statedot(4:6) = - mu * r0 / norm(r0)^3;
%     statedot(10:12) = - mu * r1 / norm(r1)^3;
% end

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = (oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3));
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

function ROE = eccentric_ROE(oe_chief, oe_deputy)
    % The orbital elements are given in the following order:
    % oe_chief = [a, e, i, omega, RAAN, M]
    ROE = zeros(1,6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(3) = oe_deputy(2) - oe_chief(2);
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = wrapToPi((oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3)));
    ROE(4) = wrapToPi(oe_deputy(4) - oe_chief(4) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(2) = wrapToPi(oe_deputy(6) - oe_chief(6) + sqrt(1 - oe_chief(2)^2) * ROE(4));
end

function P = control_gain(k, N, nu, nu_ip, nu_oop)
    P = zeros(5);
    P(1,1) = cos(nu - nu_ip)^N;
    P(2,2) = cos(nu - nu_ip)^N;
    P(3,3) = cos(nu - nu_ip)^N;
    P(4,4) = cos(nu - nu_oop)^N;
    P(5,5) = cos(nu - nu_oop)^N;
    P = P ./ k;
end

function B = modified_control_matrix(a, n, e, nu, omega)
    B = zeros(5,2);
    eta = sqrt(1-e^2);

    B(1,1) = 2 * (1+e*cos(nu)) / eta;

    B(2,1) = eta * (e + cos(nu)*(2+e*cos(nu))) / (1+e*cos(nu));

    B(3,1) = (eta / e) * sin(nu) * (2+e*cos(nu)) / (1+e*cos(nu));

    B(4,2) = eta * cos(omega + nu) / (1+e*cos(nu));

    B(5,2) = eta * sin(omega + nu) / (1+e*cos(nu));

    B = B ./ (a .* n);
end

function nu = inplane_planning(e, omega, i, delta_ex, delta_ey, delta_iy)
    syms x
    f = tan(x) * (1 + e * sin(omega) / (sin(x) * (2 + e * cos(omega) * cos(x) + e * sin(omega) * sin(x)))) / (1 + e * cos(omega) / (cos(x) * (2 + e * cos(omega) * cos(x) + e * sin(omega) * sin(x)))) - (delta_ey * tan(i) + e * cos(omega) * delta_iy) / (delta_ex * tan(i) - e * sin(omega) * delta_iy);
    nu = vpasolve(f,x,[0 2 * pi]);
end

function B = control_matrix(a, n, e, nu, omega, i)
    eta = sqrt(1 - e^2);
    B = zeros(6, 2);

    B(1, 1) = 2 * (1 + e * cos(nu)) / eta;

    B(2, 1) = eta * (eta - 1) * (2 + e * cos(nu)) * sin(nu) / (1 + e * cos(nu));

    B(3, 1) = eta * ((2 + e * cos(nu)) * cos(omega + nu) + e * cos(omega)) / (1 + e * cos(nu));
    B(3, 2) = eta * e * sin(omega) * sin(omega + nu) / (tan(i) * (1 + e * cos(nu)));

    B(4, 1) = eta * ((2 + e * cos(nu)) * sin(omega + nu) + e * sin(omega)) / (1 + e * cos(nu));
    B(4, 2) = - eta *e *cos(omega) * sin(omega + nu) / (tan(i) * (1 + e * cos(nu)));

    B(5, 2) = eta * cos(omega + nu) / (1 + e * cos(nu));

    B(6, 2) = eta * sin(omega + nu) / (1 + e * cos(nu));

    B = B ./ (a .* n);
end

function statedot = FODE_2sats(t, state, mu, control)
    statedot = zeros(size(state));
    statedot(1:3) = state(4:6);
    statedot(7:9) = state(10:12);

    r0 = state(1:3);
    r1 = state(7:9);

    statedot(4:6) = - mu * r0 / norm(r0)^3;
    statedot(10:12) = - mu * r1 / norm(r1)^3 + control;
end