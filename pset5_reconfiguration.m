%% For reconfiguration

%% Problem 1

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
R_E = 6378; % km
J2 = 0.108263e-2;

% Defining the ROE for Perigee Pass Mode (PPM)

% Initial conditions for the chief

a_chief = 36943; % km
e_chief = 0.8111;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
T = 2 * pi / n_chief;
% M_chief = n_chief * (14 * 3600 + 49 * 60);
M_chief = n_chief * (4 * 3600 + 49 * 60);
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

init_ROE_PPM = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]';
init_mod_ROE_PPM = OE2ROE_modified(oe_chief_PPM, oe_deputy_PPM)';

% Checking initial separation and initial conditions for FODE
[pos_chief, vel_chief] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief_PPM, mu);
[pos_deputy, vel_deputy] = OE2ECI(a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, nu_deputy_PPM, mu);
rho_ECI = norm(pos_deputy - pos_chief);
norm(pos_chief);

% Numerical integration without J2 perturbations
init_FODE = [pos_chief, vel_chief, pos_deputy, vel_deputy];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100);
N = 100000;
tspan = linspace(0, 100 * T, N);
[t, y] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, init_FODE, options);

mod_ROE_gt = zeros(N,6);
mod_ROE_stm = zeros(N,6);
rel_motion_RTN_gt = zeros(N,6);
rel_motion_RTN_stm = zeros(N,6);

for j=1:N
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    mod_ROE = OE2ROE_modified(oe_chief, oe_deputy);
    mod_ROE_gt(j,:) = mod_ROE;

    rel_state_ECI = y(j,7:12) - y(j,1:6);
    rel_state_RTN = ECI2RTN(y(j,1:6), rel_state_ECI, mu);
    rel_motion_RTN_gt(j, :) = rel_state_RTN;
    
    STM = STM_ROE(tspan(j), n_chief);
    mod_ROE = STM * init_mod_ROE_PPM;
    oe_deputy = ROE2OE_modified(oe_chief, mod_ROE);
    nu_chief_cur = mean2true(oe_chief(6), oe_chief(2));
    nu_deputy_cur = mean2true(oe_deputy(6), oe_deputy(2));
    [pos_c, vel_c] = OE2ECI(oe_chief(1), oe_chief(2), oe_chief(3), oe_chief(4), oe_chief(5), nu_chief_cur, mu);
    [pos_d, vel_d] = OE2ECI(oe_deputy(1), oe_deputy(2), oe_deputy(3), oe_deputy(4), oe_deputy(5), nu_deputy_cur, mu);
    rel_state_ECI = [pos_d', vel_d'] - [pos_c', vel_c'];
    rel_state_RTN = ECI2RTN([pos_c', vel_c'], rel_state_ECI, mu);
    rel_motion_RTN_stm(j, :) = rel_state_RTN;
    mod_ROE_stm(j, :) = mod_ROE;
end

% figure
% subplot(3,2,1)
% hold on
% plot(t / T, mod_ROE_gt(:, 1) * a_chief * 1e3)
% plot(t / T, mod_ROE_stm(:, 1) * a_chief * 1e3)
% hold off
% grid on
% legend('Ground Truth', 'STM')
% ylabel('\deltaa [m]')
% 
% subplot(3,2,3)
% hold on
% plot(t / T, mod_ROE_gt(:, 2) * a_chief * 1e3)
% plot(t / T, mod_ROE_stm(:, 2) * a_chief * 1e3)
% hold off
% grid on
% ylabel('\delta\lambda [m]')
% 
% subplot(3,2,5)
% hold on
% plot(t / T, mod_ROE_gt(:, 3) * a_chief * 1e3)
% plot(t / T, mod_ROE_stm(:, 3) * a_chief * 1e3)
% hold off
% grid on
% ylabel("\deltae_x' [m]")
% xlabel('Orbital Periods')
% 
% subplot(3,2,2)
% hold on
% plot(t / T, mod_ROE_gt(:, 4) * a_chief * 1e3)
% plot(t / T, mod_ROE_stm(:, 4) * a_chief * 1e3)
% hold off
% grid on
% ylabel("\deltae_y' [m]")
% 
% subplot(3,2,4)
% hold on
% plot(t / T, mod_ROE_gt(:, 5) * a_chief * 1e3)
% plot(t / T, mod_ROE_stm(:, 5) * a_chief * 1e3)
% hold off
% grid on
% ylabel('\deltai_x [m]')
% 
% subplot(3,2,6)
% hold on
% plot(t / T, mod_ROE_gt(:, 6) * a_chief * 1e3)
% plot(t / T, mod_ROE_stm(:, 6) * a_chief * 1e3)
% hold off
% grid on
% ylabel('\deltai_y [m]')
% xlabel('Orbital Periods')
% 
% figure
% subplot(3,1,1)
% hold on
% plot(mod_ROE_gt(:, 2) * a_chief * 1e3, mod_ROE_gt(:, 1) * a_chief * 1e3)
% plot(mod_ROE_stm(:, 2) * a_chief * 1e3, mod_ROE_stm(:, 1) * a_chief * 1e3, 'o')
% hold off
% axis equal
% grid on
% legend('Ground Truth', 'STM')
% xlabel('\delta\lambda [m]')
% ylabel('\deltaa [m]')
% 
% subplot(3,1,2)
% hold on
% plot(mod_ROE_gt(:,3) * a_chief * 1e3, mod_ROE_gt(:,4) * a_chief * 1e3)
% plot(mod_ROE_stm(:,3) * a_chief * 1e3, mod_ROE_stm(:,4) * a_chief * 1e3, 'o')
% hold off
% axis equal
% grid on
% xlabel("\deltae_x' [m]")
% ylabel("\deltae_y' [m]")
% 
% subplot(3,1,3)
% hold on
% plot(mod_ROE_gt(:,5) * a_chief * 1e3, mod_ROE_gt(:,6) * a_chief * 1e3)
% plot(mod_ROE_stm(:,5) * a_chief * 1e3, mod_ROE_stm(:,6) * a_chief * 1e3,'o')
% hold off
% axis equal
% grid on
% xlabel('\deltai_x [m]')
% ylabel('\deltai_y [m]')
% 
% figure
% subplot(2,2,1)
% plot(rel_motion_RTN_gt(:,2), rel_motion_RTN_gt(:,1))
% axis equal
% grid on
% xlabel('T-axis [km]')
% ylabel('R-axis [km]')
% 
% subplot(2,2,2)
% plot(rel_motion_RTN_gt(:,3), rel_motion_RTN_gt(:,1))
% axis equal
% grid on
% xlabel('N-axis [km]')
% ylabel('R-axis [km]')
% 
% subplot(2,2,3)
% plot(rel_motion_RTN_gt(:,2), rel_motion_RTN_gt(:,3))
% axis equal
% grid on
% xlabel('T-axis [km]')
% ylabel('N-axis [km]')
% 
% subplot(2,2,4)
% plot3(rel_motion_RTN_gt(:,1), rel_motion_RTN_gt(:,2), rel_motion_RTN_gt(:,3))
% view(3)
% axis equal
% grid on
% xlabel('R-axis [km]')
% ylabel('T-axis [km]')
% zlabel('N-axis [km]')
% 
% figure
% subplot(2,2,1)
% plot(rel_motion_RTN_gt(:,5), rel_motion_RTN_gt(:,4))
% axis equal
% grid on
% xlabel('T-axis [km/s]')
% ylabel('R-axis [km/s]')
% 
% subplot(2,2,2)
% plot(rel_motion_RTN_gt(:,6), rel_motion_RTN_gt(:,4))
% axis equal
% grid on
% xlabel('N-axis [km/s]')
% ylabel('R-axis [km/s]')
% 
% subplot(2,2,3)
% plot(rel_motion_RTN_gt(:,5), rel_motion_RTN_gt(:,6))
% axis equal
% grid on
% xlabel('T-axis [km/s]')
% ylabel('N-axis [km/s]')
% 
% subplot(2,2,4)
% plot3(rel_motion_RTN_gt(:,4), rel_motion_RTN_gt(:,5), rel_motion_RTN_gt(:,6))
% view(3)
% axis equal
% grid on
% xlabel('R-axis [km/s]')
% ylabel('T-axis [km/s]')
% zlabel('N-axis [km/S]')
% 
% figure
% subplot(3,2,1)
% plot(t / T, abs(rel_motion_RTN_gt(:,1) - rel_motion_RTN_stm(:,1)))
% ylabel('Eror in R-axis [km]')
% 
% subplot(3,2,3)
% plot(t / T, abs(rel_motion_RTN_gt(:,2) - rel_motion_RTN_stm(:,2)))
% ylabel('Error in T-axis [km]')
% 
% subplot(3,2,5)
% plot(t / T, abs(rel_motion_RTN_gt(:,3) - rel_motion_RTN_stm(:,3)))
% ylabel('Error in N-axis [km]')
% xlabel("Orbital Periods")
% 
% subplot(3,2,2)
% plot(t / T, abs(rel_motion_RTN_gt(:,4) - rel_motion_RTN_stm(:,4)))
% ylabel('Error in R-axis [km/s]')
% 
% subplot(3,2,4)
% plot(t / T, abs(rel_motion_RTN_gt(:,5) - rel_motion_RTN_stm(:,5)))
% ylabel('Error in T-axis [km/s]')
% 
% subplot(3,2,6)
% plot(t / T, abs(rel_motion_RTN_gt(:,6) - rel_motion_RTN_stm(:,6)))
% ylabel('Error in N-axis [km/s]')
% xlabel('Orbital Periods')

%% Problem 2
% Reconfiguration from Perigee Pass Mode to Inertial Attitude Mode using least squares

% Defining the ROE for Inertial Attitude Mode (IAM)

% Initial conditions for the chief
M_chief = n_chief * (6 * 3600 + 49 * 60);
nu_chief_IAM = mean2true(M_chief, e_chief);
oe_chief_IAM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

% Initial conditions for the deputy

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
oe_deputy_IAM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

init_mod_ROE_IAM = OE2ROE_modified(oe_chief_IAM, oe_deputy_IAM)';

% Compute the pseudostate
STM_f0 = STM_ROE(2 * 3600, n_chief); % 2 hours to perform the reconfiguration
delta_alpha = init_mod_ROE_IAM - STM_f0 * init_mod_ROE_PPM;

% Compuuting delta-V lower bound
abs(delta_alpha(1)) * n_chief * a_chief * sqrt(1-e_chief^2) / (2*(1+e_chief))
abs(delta_alpha(2)) * n_chief * a_chief * sqrt(1-e_chief^2) / (3*(1+e_chief)*n_chief*2*3600)
norm([delta_ex - init_ROE_PPM(3), delta_ey - init_ROE_PPM(4)]) * n_chief * a_chief * sqrt(1-e_chief^2) / sqrt(3*e_chief^4 - 7 * e_chief^2 + 4)

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

% super_STM = [STM_f2(:,1), STM_f4(:,1), STM_f2(:,2), STM_f4(:,2), STM_f2(:,3), STM_f4(:,3), STM_f2(:,4), STM_f4(:,4), STM_f2(:,5), STM_f4(:,5), STM_f2(:,6), STM_f4(:,6)];
% super_control = [STM_f0*B_0, STM_f1*B_1, STM_f2*B_2, STM_f3*B_3, STM_f4*B_4, STM_f5*B_5, STM_f6*B_6];
super_control = [STM_f0*B_0, STM_f1*B_1, STM_f2*B_2, STM_f3*B_3, STM_f4*B_4, STM_f5*B_5];
delta_v = (super_control' * super_control)^(-1) * super_control' * delta_alpha; % should be in km/s

% Applying the maneuvers

tspan_01 = linspace(4*3600 + 49 * 60, 5 * 3600 + 9 * 60, 1000);
[pos_chief_0, vel_chief_0] = OE2ECI(oe_chief_PPM(1), oe_chief_PPM(2), oe_chief_PPM(3), oe_chief_PPM(4), oe_chief_PPM(5), nu_chief_PPM, mu);
[pos_deputy_0, vel_deputy_0] = OE2ECI(oe_deputy_PPM(1), oe_deputy_PPM(2), oe_deputy_PPM(3), oe_deputy_PPM(4), oe_deputy_PPM(5), nu_deputy_PPM, mu);
% RTN unit vectors
R_0 = pos_chief' / norm(pos_chief);
h_0 = cross(pos_chief, vel_chief)';
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
    ROE = OE2ROE_modified(oe_chief, oe_deputy);
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
    ROE = OE2ROE_modified(oe_chief, oe_deputy);
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
    ROE = OE2ROE_modified(oe_chief, oe_deputy);
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
    ROE = OE2ROE_modified(oe_chief, oe_deputy);
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
    ROE = OE2ROE_modified(oe_chief, oe_deputy);
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
    ROE = OE2ROE_modified(oe_chief, oe_deputy);
    ROE_gt_56(j,:) = ROE;
end

figure
subplot(3,2,1)
hold on
plot((4 * 3600 + 49 * 60) / T, init_mod_ROE_PPM(1) * a_chief * 1e3, 'ro')
plot(t_56(end) / T, init_mod_ROE_IAM(1) * a_chief * 1e3, 'rx')
plot(t_01 / T, ROE_gt_01(:, 1) * a_chief * 1e3, 'b-')
plot(t_12 / T, ROE_gt_12(:, 1) * a_chief * 1e3, 'b-')
plot(t_23 / T, ROE_gt_23(:, 1) * a_chief * 1e3, 'b-')
plot(t_34 / T, ROE_gt_34(:, 1) * a_chief * 1e3, 'b-')
plot(t_45 / T, ROE_gt_45(:, 1) * a_chief * 1e3, 'b-')
plot(t_56 / T, ROE_gt_56(:, 1) * a_chief * 1e3, 'b-')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
legend('Initial', 'Target')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot((4 * 3600 + 49 * 60) / T, init_mod_ROE_PPM(2) * a_chief * 1e3, 'ro')
plot(t_01 / T, ROE_gt_01(:, 2) * a_chief * 1e3, 'b-')
plot(t_12 / T, ROE_gt_12(:, 2) * a_chief * 1e3, 'b-')
plot(t_23 / T, ROE_gt_23(:, 2) * a_chief * 1e3, 'b-')
plot(t_34 / T, ROE_gt_34(:, 2) * a_chief * 1e3, 'b-')
plot(t_45 / T, ROE_gt_45(:, 2) * a_chief * 1e3, 'b-')
plot(t_56 / T, ROE_gt_56(:, 2) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, init_mod_ROE_IAM(2) * a_chief * 1e3, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot((4 * 3600 + 49 * 60) / T, init_mod_ROE_PPM(3) * a_chief * 1e3, 'ro')
plot(t_01 / T, ROE_gt_01(:, 3) * a_chief * 1e3, 'b-')
plot(t_12 / T, ROE_gt_12(:, 3) * a_chief * 1e3, 'b-')
plot(t_23 / T, ROE_gt_23(:, 3) * a_chief * 1e3, 'b-')
plot(t_34 / T, ROE_gt_34(:, 3) * a_chief * 1e3, 'b-')
plot(t_45 / T, ROE_gt_45(:, 3) * a_chief * 1e3, 'b-')
plot(t_56 / T, ROE_gt_56(:, 3) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, init_mod_ROE_IAM(3) * a_chief * 1e3, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k')
hold off
grid on
ylabel("\deltae_x' [m]")
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot((4 * 3600 + 49 * 60) / T, init_mod_ROE_PPM(4) * a_chief * 1e3, 'ro')
plot(t_01 / T, ROE_gt_01(:, 4) * a_chief * 1e3, 'b-')
plot(t_12 / T, ROE_gt_12(:, 4) * a_chief * 1e3, 'b-')
plot(t_23 / T, ROE_gt_23(:, 4) * a_chief * 1e3, 'b-')
plot(t_34 / T, ROE_gt_34(:, 4) * a_chief * 1e3, 'b-')
plot(t_45 / T, ROE_gt_45(:, 4) * a_chief * 1e3, 'b-')
plot(t_56 / T, ROE_gt_56(:, 4) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, init_mod_ROE_IAM(4) * a_chief * 1e3, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
ylabel("\deltae_y' [m]")

subplot(3,2,4)
hold on
plot((4 * 3600 + 49 * 60) / T, init_mod_ROE_PPM(5) * a_chief * 1e3, 'ro')
plot(t_01 / T, ROE_gt_01(:, 5) * a_chief * 1e3, 'b-')
plot(t_12 / T, ROE_gt_12(:, 5) * a_chief * 1e3, 'b-')
plot(t_23 / T, ROE_gt_23(:, 5) * a_chief * 1e3, 'b-')
plot(t_34 / T, ROE_gt_34(:, 5) * a_chief * 1e3, 'b-')
plot(t_45 / T, ROE_gt_45(:, 5) * a_chief * 1e3, 'b-')
plot(t_56 / T, ROE_gt_56(:, 5) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, init_mod_ROE_IAM(5) * a_chief * 1e3, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot((4 * 3600 + 49 * 60) / T, init_mod_ROE_PPM(6) * a_chief * 1e3, 'ro')
plot(t_01 / T, ROE_gt_01(:, 6) * a_chief * 1e3, 'b-')
plot(t_12 / T, ROE_gt_12(:, 6) * a_chief * 1e3, 'b-')
plot(t_23 / T, ROE_gt_23(:, 6) * a_chief * 1e3, 'b-')
plot(t_34 / T, ROE_gt_34(:, 6) * a_chief * 1e3, 'b-')
plot(t_45 / T, ROE_gt_45(:, 6) * a_chief * 1e3, 'b-')
plot(t_56 / T, ROE_gt_56(:, 6) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, init_mod_ROE_IAM(6) * a_chief * 1e3, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k')
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

% Control errors and maneuver planning

figure
subplot(3,2,1)
hold on
plot((4 * 3600 + 49 * 60) / T, (init_mod_ROE_PPM(1)-init_mod_ROE_IAM(1)) * a_chief * 1e3, 'ro')
plot(t_56(end) / T, 0, 'rx')
plot(t_01 / T, (ROE_gt_01(:, 1)-init_mod_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_12 / T, (ROE_gt_12(:, 1)-init_mod_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_23 / T, (ROE_gt_23(:, 1)-init_mod_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_34 / T, (ROE_gt_34(:, 1)-init_mod_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_45 / T, (ROE_gt_45(:, 1)-init_mod_ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_56 / T, (ROE_gt_56(:, 1)-init_mod_ROE_IAM(1)) * a_chief * 1e3, 'b-')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
legend('Initial', 'Target')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot((4 * 3600 + 49 * 60) / T, (init_mod_ROE_PPM(2)-init_mod_ROE_IAM(2)) * a_chief * 1e3, 'ro')
plot(t_01 / T, (ROE_gt_01(:, 2)-init_mod_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_12 / T, (ROE_gt_12(:, 2)-init_mod_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_23 / T, (ROE_gt_23(:, 2)-init_mod_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_34 / T, (ROE_gt_34(:, 2)-init_mod_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_45 / T, (ROE_gt_45(:, 2)-init_mod_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_56 / T, (ROE_gt_56(:, 2)-init_mod_ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, 0, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot((4 * 3600 + 49 * 60) / T, (init_mod_ROE_PPM(3)-init_mod_ROE_IAM(3)) * a_chief * 1e3, 'ro')
plot(t_01 / T, (ROE_gt_01(:, 3)-init_mod_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_12 / T, (ROE_gt_12(:, 3)-init_mod_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_23 / T, (ROE_gt_23(:, 3)-init_mod_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_34 / T, (ROE_gt_34(:, 3)-init_mod_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_45 / T, (ROE_gt_45(:, 3)-init_mod_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_56 / T, (ROE_gt_56(:, 3)-init_mod_ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, 0, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k')
hold off
grid on
ylabel("\deltae_x' [m]")
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot((4 * 3600 + 49 * 60) / T, (init_mod_ROE_PPM(4)-init_mod_ROE_IAM(4)) * a_chief * 1e3, 'ro')
plot(t_01 / T, (ROE_gt_01(:, 4)-init_mod_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_12 / T, (ROE_gt_12(:, 4)-init_mod_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_23 / T, (ROE_gt_23(:, 4)-init_mod_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_34 / T, (ROE_gt_34(:, 4)-init_mod_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_45 / T, (ROE_gt_45(:, 4)-init_mod_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_56 / T, (ROE_gt_56(:, 4)-init_mod_ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, 0, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
ylabel("\deltae_y' [m]")

subplot(3,2,4)
hold on
plot((4 * 3600 + 49 * 60) / T, (init_mod_ROE_PPM(5)-init_mod_ROE_IAM(5)) * a_chief * 1e3, 'ro')
plot(t_01 / T, (ROE_gt_01(:, 5)-init_mod_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_12 / T, (ROE_gt_12(:, 5)-init_mod_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_23 / T, (ROE_gt_23(:, 5)-init_mod_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_34 / T, (ROE_gt_34(:, 5)-init_mod_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_45 / T, (ROE_gt_45(:, 5)-init_mod_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_56 / T, (ROE_gt_56(:, 5)-init_mod_ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, 0, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot((4 * 3600 + 49 * 60) / T, (init_mod_ROE_PPM(6)-init_mod_ROE_IAM(6)) * a_chief * 1e3, 'ro')
plot(t_01 / T, (ROE_gt_01(:, 6)-init_mod_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_12 / T, (ROE_gt_12(:, 6)-init_mod_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_23 / T, (ROE_gt_23(:, 6)-init_mod_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_34 / T, (ROE_gt_34(:, 6)-init_mod_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_45 / T, (ROE_gt_45(:, 6)-init_mod_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_56 / T, (ROE_gt_56(:, 6)-init_mod_ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_56(end) / T, 0, 'rx')
xline([t_01(1)/T t_01(end)/T t_12(end)/T t_23(end)/T t_34(end)/T t_45(end)/T], '--k')
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

%% Funtions

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = wrapToPi((oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3)));
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

function oe_deputy = ROE2OE_modified(oe_chief, modified_ROE)
    oe_deputy = zeros(1,6);
    oe_deputy(1) = oe_chief(1) + oe_chief(1) * modified_ROE(1);
    oe_deputy(2) = oe_chief(2) + modified_ROE(3);
    oe_deputy(6) = modified_ROE(2) - modified_ROE(4) + oe_chief(6);
    oe_deputy(3) = oe_chief(3) + modified_ROE(5);
    oe_deputy(5) = oe_chief(5) + modified_ROE(6) / sin(oe_chief(3));
    oe_deputy(4) = modified_ROE(4) + oe_chief(4) - (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3));
end

function STM = STM_ROE(t, n)
    STM = eye(6);
    STM(2,1) = -1.5 * n * t;
end

function B = control_matrix(nu, e, a, n, omega)
    eta = sqrt(1-e^2);                   
    B_nd = zeros(6, 3);
    B_nd(1,1) = 2 * e * sin(nu) / eta;
    B_nd(1,2) = (1+e*cos(nu)) * (2 / eta);

    B_nd(2,1) = - 2 * eta^2 / (1+e*cos(nu));
%     B_nd(2,1) = eta * (eta - 1) * cos(nu) / e - (2 * eta^2) / (1+e*cos(nu));

    B_nd(3,1) = eta * sin(nu);
%     B(3,2) = ((2+e*cos(nu))*cos(nu) + e) / (1+e*cos(nu)) * eta;
    B_nd(3,2) = (eta / (1+e*cos(nu))) * (e + cos(nu) * (2+e*cos(nu))); % these
%     two lines different somehow

    B_nd(4,1) = - eta * cos(nu) / e;
%     B(4,2) = ((2+e*cos(nu))*sin(nu)) / (1+e*cos(nu)) * eta / e;
    B_nd(4,2) = (eta / e) * sin(nu) * (2+e*cos(nu)) / (1+e*cos(nu));

    B_nd(5,3) = eta * cos(nu+omega) / (1+e*cos(nu));
    B_nd(6,3) = eta * sin(nu+omega) / (1+e*cos(nu));
    B = B_nd / (n * a);
end

function B = control_matrix_2(nu, e, i, a, n, omega)
    eta = sqrt(1-e^2);                   
    B_nd = zeros(6, 3);
    B_nd(1,1) = 2 * e * sin(nu) / eta;
    B_nd(1,2) = (1+e*cos(nu)) * 2 / eta;
    B_nd(2,1) = - 2 * eta^2 / (1+e*cos(nu));
    B_nd(3,1) = eta * sin(nu + omega);
    B_nd(3,2) = ((2+e*cos(nu))*cos(nu+omega) + e*cos(omega)) / (1+e*cos(nu)) * eta;
    B_nd(3,3) = eta * e * sin(omega) * sin(nu+omega) / ((tan(i)) * (1+e*cos(nu)));
    B_nd(4,1) = - eta * cos(nu+omega);
    B_nd(4,2) = ((2+e*cos(nu))*sin(nu+omega) + e*sin(omega)) / (1+e*cos(nu)) * eta;
    B_nd(4,3) = - eta * e * cos(omega) * sin(nu+omega) / ((tan(i)) * (1+e*cos(nu)));
    B_nd(5,3) = eta * cos(nu+omega) / (1+e*cos(nu));
    B_nd(6,3) = eta * sin(nu+omega) / (1+e*cos(nu));
    B = B_nd / (n * a);
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