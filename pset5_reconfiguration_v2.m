close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100); % options for numerical integration

% Computing the initial modified ROE for the reconfiguration (Perigee Pass
% Mode)

a_chief = 36943; % km
e_chief = 0.8111;
% e_chief = 1e-5;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
T_chief = 2 * pi / n_chief;
M_chief = n_chief * (4 * 3600 + 49 * 60); % reconfiguration starts at t=4h49
nu_chief_PPM = mean2true(M_chief, e_chief);
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
nu_deputy_PPM = mean2true(M_deputy, e_deputy);
oe_deputy_PPM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

modified_ROE_PPM = eccentric_ROE(oe_chief_PPM, oe_deputy_PPM);
ROE_PPM = OE2ROE(oe_chief_PPM, oe_deputy_PPM);

% Computing the final modified ROE for the reconfiguration (Inertial
% Attitude Mode)

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

% Compute the pseudostate
STM_f0 = STM_ROE(2 * 3600, n_chief); % 2 hours to perform the reconfiguration
delta_alpha = modified_ROE_IAM' - STM_f0 * modified_ROE_PPM'; % column vector

% % Compuuting delta-V lower bound
% abs(delta_alpha(1)) * n_chief * a_chief * sqrt(1-e_chief^2) / (2*(1+e_chief))
% abs(delta_alpha(2)) * n_chief * a_chief * sqrt(1-e_chief^2) / (3*(1+e_chief)*n_chief*2*3600)
% norm([delta_ex - ROE_PPM(3), delta_ey - ROE_PPM(4)]) * n_chief * a_chief * sqrt(1-e_chief^2) / sqrt(3*e_chief^4 - 7 * e_chief^2 + 4)

% Seeting the least squares to solve the control problem
B_0 = control_matrix(mean2true(oe_chief_IAM(6), e_chief), e_chief, a_chief, n_chief, omega_chief);
STM_f1 = STM_ROE(3600 + 40 * 60, n_chief);
B_1 = control_matrix(mean2true((5*3600+9*60)*n_chief, e_chief), e_chief, a_chief, n_chief, omega_chief);
STM_f2 = STM_ROE(3600 + 20 * 60, n_chief);
B_2 = control_matrix(mean2true((5*3600+29*60)*n_chief, e_chief), e_chief, a_chief, n_chief, omega_chief);
STM_f3 = STM_ROE(3600, n_chief);
B_3 = control_matrix(mean2true((5*3600+49*60)*n_chief, e_chief), e_chief, a_chief, n_chief, omega_chief);
STM_f4 = STM_ROE(40 * 60, n_chief);
B_4 = control_matrix(mean2true((6*3600+9*60)*n_chief, e_chief), e_chief, a_chief, n_chief, omega_chief);
STM_f5 = STM_ROE(20 * 60, n_chief);
B_5 = control_matrix(mean2true((6*3600+29*60)*n_chief, e_chief), e_chief, a_chief, n_chief, omega_chief);
STM_f6 = eye(6);
B_6 = control_matrix(mean2true((6*3600+49*60)*n_chief, e_chief), e_chief, a_chief, n_chief, omega_chief);

super_control = [STM_f0*B_0, STM_f1*B_1, STM_f2*B_2, STM_f3*B_3, STM_f4*B_4, STM_f5*B_5];
delta_v = (super_control' * super_control)^(-1) * super_control' * delta_alpha; % should be in km/s

% Applying the maneuvers

tspan_01 = linspace(4*3600 + 49 * 60, 5 * 3600 + 9 * 60, 1000);
[pos_chief_0, vel_chief_0] = OE2ECI(oe_chief_PPM(1), oe_chief_PPM(2), oe_chief_PPM(3), oe_chief_PPM(4), oe_chief_PPM(5), nu_chief_PPM, mu);
[pos_deputy_0, vel_deputy_0] = OE2ECI(oe_deputy_PPM(1), oe_deputy_PPM(2), oe_deputy_PPM(3), oe_deputy_PPM(4), oe_deputy_PPM(5), nu_deputy_PPM, mu);

% RTN unit vectors
R_0 = pos_chief_0' / norm(pos_chief_0);
h_0 = cross(pos_chief_0, vel_chief_0)';
N_0 = h_0 / norm(h_0);
T_0 = cross(N_0, R_0);

% Rotation matrix from ECI to RTN
rotation_0 = [R_0; T_0; N_0];

vel_deputy_0 = vel_deputy_0' + (rotation_0'*delta_v(1:3))';
init_cond_01 = [pos_chief_0', vel_chief_0', pos_deputy_0', vel_deputy_0];

[t_01, y_01] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_01, init_cond_01, options);

ROE_gt_01 = zeros(1000,6);

for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_01(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_01(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = eccentric_ROE(oe_chief, oe_deputy);
    ROE_gt_01(j,:) = ROE;
end

% RTN unit vectors
R_1 = y_01(end, 1:3) / norm(y_01(end, 1:3));
h_1 = cross(y_01(end, 1:3), y_01(end, 4:6));
N_1 = h_1 / norm(h_1);
T_1 = cross(N_1, R_1);

% Rotation matrix from ECI to RTN
rotation_1 = [R_1; T_1; N_1];

vel_deputy_1 = y_01(end, 10:12) + (rotation_1'*delta_v(4:6))';
init_cond_12 = [y_01(end, 1:9) vel_deputy_1];
tspan_12 = linspace(5 * 3600 + 9 * 60, 5 * 3600 + 29 * 60, 1000);
[t_12, y_12] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_12, init_cond_12, options);

ROE_gt_12 = zeros(1000,6);

for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_12(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_12(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = eccentric_ROE(oe_chief, oe_deputy);
    ROE_gt_12(j,:) = ROE;
end

% RTN unit vectors
R_2 = y_12(end, 1:3) / norm(y_12(end, 1:3));
h_2 = cross(y_12(end, 1:3), y_12(end, 4:6));
N_2 = h_2 / norm(h_2);
T_2 = cross(N_2, R_2);

% Rotation matrix from ECI to RTN
rotation_2 = [R_2; T_2; N_2];

vel_deputy_2 = y_12(end, 10:12) + (rotation_2'*delta_v(7:9))';
init_cond_23 = [y_12(end, 1:9) vel_deputy_2];
tspan_23 = linspace(5 * 3600 + 29 * 60, 5 * 3600 + 49 * 60, 1000);
[t_23, y_23] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_23, init_cond_23, options);

ROE_gt_23 = zeros(1000,6);

for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_23(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_23(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = eccentric_ROE(oe_chief, oe_deputy);
    ROE_gt_23(j,:) = ROE;
end

% RTN unit vectors
R_3 = y_23(end, 1:3) / norm(y_23(end, 1:3));
h_3 = cross(y_23(end, 1:3), y_23(end, 4:6));
N_3 = h_3 / norm(h_3);
T_3 = cross(N_3, R_3);

% Rotation matrix from ECI to RTN
rotation_3 = [R_3; T_3; N_3];

vel_deputy_3 = y_23(end, 10:12) + (rotation_3'*delta_v(10:12))';
init_cond_34 = [y_23(end, 1:9) vel_deputy_3];
tspan_34 = linspace(5 * 3600 + 49 * 60, 6 * 3600 + 9 * 60, 1000);
[t_34, y_34] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_34, init_cond_34, options);

ROE_gt_34 = zeros(1000,6);

for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_34(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_34(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = eccentric_ROE(oe_chief, oe_deputy);
    ROE_gt_34(j,:) = ROE;
end

% RTN unit vectors
R_4 = y_34(end, 1:3) / norm(y_34(end, 1:3));
h_4 = cross(y_34(end, 1:3), y_34(end, 4:6));
N_4 = h_4 / norm(h_4);
T_4 = cross(N_4, R_4);

% Rotation matrix from ECI to RTN
rotation_4 = [R_4; T_4; N_4];

vel_deputy_4 = y_34(end, 10:12) + (rotation_4'*delta_v(13:15))';
init_cond_45 = [y_34(end, 1:9) vel_deputy_4];
tspan_45 = linspace(6 * 3600 + 9 * 60, 6 * 3600 + 29 * 60, 1000);
[t_45, y_45] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_45, init_cond_45, options);

ROE_gt_45 = zeros(1000,6);

for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_45(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_45(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = eccentric_ROE(oe_chief, oe_deputy);
    ROE_gt_45(j,:) = ROE;
end

% RTN unit vectors
R_5 = y_45(end, 1:3) / norm(y_45(end, 1:3));
h_5 = cross(y_45(end, 1:3), y_45(end, 4:6));
N_5 = h_5 / norm(h_5);
T_5 = cross(N_5, R_5);

% Rotation matrix from ECI to RTN
rotation_5 = [R_5; T_5; N_5];

vel_deputy_5 = y_45(end, 10:12) + (rotation_5'*delta_v(16:18))';
init_cond_56 = [y_45(end, 1:9) vel_deputy_5];
tspan_56 = linspace(6 * 3600 + 29 * 60, 6 * 3600 + 49 * 60, 1000);
[t_56, y_56] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_56, init_cond_56, options);

ROE_gt_56 = zeros(1000,6);

for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_56(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_56(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = eccentric_ROE(oe_chief, oe_deputy);
    ROE_gt_56(j,:) = ROE;
end

figure
subplot(3,2,1)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, modified_ROE_PPM(1) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, modified_ROE_IAM(1) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, ROE_gt_01(:, 1) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_gt_12(:, 1) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_gt_23(:, 1) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_gt_34(:, 1) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_gt_45(:, 1) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_gt_56(:, 1) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
legend('Initial', 'Target')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, modified_ROE_PPM(2) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_gt_01(:, 2) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_gt_12(:, 2) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_gt_23(:, 2) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_gt_34(:, 2) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_gt_45(:, 2) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_gt_56(:, 2) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, modified_ROE_IAM(2) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, modified_ROE_PPM(3) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_gt_01(:, 3) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_gt_12(:, 3) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_gt_23(:, 3) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_gt_34(:, 3) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_gt_45(:, 3) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_gt_56(:, 3) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, modified_ROE_IAM(3) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel("\deltae_x' [m]")
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, modified_ROE_PPM(4) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_gt_01(:, 4) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_gt_12(:, 4) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_gt_23(:, 4) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_gt_34(:, 4) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_gt_45(:, 4) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_gt_56(:, 4) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, modified_ROE_IAM(4) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
ylabel("\deltae_y' [m]")

subplot(3,2,4)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, modified_ROE_PPM(5) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_gt_01(:, 5) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_gt_12(:, 5) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_gt_23(:, 5) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_gt_34(:, 5) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_gt_45(:, 5) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_gt_56(:, 5) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, modified_ROE_IAM(5) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, modified_ROE_PPM(6) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_gt_01(:, 6) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_gt_12(:, 6) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_gt_23(:, 6) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_gt_34(:, 6) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_gt_45(:, 6) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_gt_56(:, 6) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, modified_ROE_IAM(6) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

% Control errors and maneuver planning

figure
subplot(3,2,1)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, (modified_ROE_PPM(1)-modified_ROE_IAM(1)) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, 0, 'rx')
plot(t_01 / T_chief, (ROE_gt_01(:, 1)-modified_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_gt_12(:, 1)-modified_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_gt_23(:, 1)-modified_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_gt_34(:, 1)-modified_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_gt_45(:, 1)-modified_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_gt_56(:, 1)-modified_ROE_IAM(1)) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
legend('Initial', 'Target')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, (modified_ROE_PPM(2)-modified_ROE_IAM(2)) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, (ROE_gt_01(:, 2)-modified_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_gt_12(:, 2)-modified_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_gt_23(:, 2)-modified_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_gt_34(:, 2)-modified_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_gt_45(:, 2)-modified_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_gt_56(:, 2)-modified_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, 0, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, (modified_ROE_PPM(3)-modified_ROE_IAM(3)) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, (ROE_gt_01(:, 3)-modified_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_gt_12(:, 3)-modified_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_gt_23(:, 3)-modified_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_gt_34(:, 3)-modified_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_gt_45(:, 3)-modified_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_gt_56(:, 3)-modified_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, 0, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel("\deltae_x' [m]")
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, (modified_ROE_PPM(4)-modified_ROE_IAM(4)) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, (ROE_gt_01(:, 4)-modified_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_gt_12(:, 4)-modified_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_gt_23(:, 4)-modified_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_gt_34(:, 4)-modified_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_gt_45(:, 4)-modified_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_gt_56(:, 4)-modified_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, 0, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
ylabel("\deltae_y' [m]")

subplot(3,2,4)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, (modified_ROE_PPM(5)-modified_ROE_IAM(5)) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, (ROE_gt_01(:, 5)-modified_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_gt_12(:, 5)-modified_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_gt_23(:, 5)-modified_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_gt_34(:, 5)-modified_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_gt_45(:, 5)-modified_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_gt_56(:, 5)-modified_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, 0, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, (modified_ROE_PPM(6)-modified_ROE_IAM(6)) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, (ROE_gt_01(:, 6)-modified_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_gt_12(:, 6)-modified_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_gt_23(:, 6)-modified_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_gt_34(:, 6)-modified_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_gt_45(:, 6)-modified_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_gt_56(:, 6)-modified_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, 0, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

%% Functions

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

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end

function STM = STM_ROE(t, n)
    STM = eye(6);
    STM(2,1) = -1.5 * n * t;
end

function B = control_matrix(nu,e, a, n, omega)
    B = zeros(6,3);
    eta = sqrt(1-e^2);

    B(1,1) = 2 / eta * e * sin(nu);
    B(1,2) = 2 * (1+e*cos(nu)) / eta;

    B(2,1) = -2 * eta^2 / (1+e*cos(nu));

    B(3,1) = eta * sin(nu);
    B(3,2) = eta * (e + cos(nu)*(2+e*cos(nu))) / (1+e*cos(nu));

    B(4,1) = - eta / e * cos(nu);
    B(4,2) = (eta / e) * sin(nu) * (2+e*cos(nu)) / (1+e*cos(nu));

    B(5,3) = eta * cos(omega + nu) / (1+e*cos(nu));

    B(6,3) = eta * sin(omega + nu) / (1+e*cos(nu));

    B = B ./ (a .* n);
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