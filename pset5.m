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
M_chief = n_chief * (14 * 3600 + 49 * 60);
nu_chief = mean2true(M_chief, e_chief);
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
nu_deputy = mean2true(M_deputy, e_deputy);
oe_deputy_PPM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

init_mod_ROE_PPM = OE2ROE_modified(oe_chief_PPM, oe_deputy_PPM)';

% Checking initial separation and initial conditions for FODE
[pos_chief, vel_chief] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief, mu);
[pos_deputy, vel_deputy] = OE2ECI(a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, nu_deputy, mu);
rho_ECI = norm(pos_deputy - pos_chief)
norm(pos_chief);

% Numerical integration without J2 perturbations
init_FODE = [pos_chief, vel_chief, pos_deputy, vel_deputy];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100);
N = 100000;
tspan = linspace(0, 50 * T, N);
[t, y] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, init_FODE, options);

mod_ROE_gt = zeros(N,6);
mod_ROE_stm = zeros(N,6);
rel_motion_RTN = zeros(N,6);

for j=1:N
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    mod_ROE = OE2ROE_modified(oe_chief, oe_deputy);
    mod_ROE_gt(j,:) = mod_ROE;

    rel_state_ECI = y(j,7:12) - y(j,1:6);
    rel_state_RTN = ECI2RTN(y(j,1:6), rel_state_ECI, mu);
    rel_motion_RTN(j, :) = rel_state_RTN;
    
    STM = STM_ROE(tspan(j), n_chief);
    mod_ROE = STM * init_mod_ROE_PPM;
    mod_ROE_stm(j, :) = mod_ROE;
end

figure
subplot(3,2,1)
hold on
plot(t / T, mod_ROE_gt(:, 1) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 1) * a_chief * 1e3)
hold off
grid on
legend('Ground Truth', 'STM')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(t / T, mod_ROE_gt(:, 2) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 2) * a_chief * 1e3)
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(t / T, mod_ROE_gt(:, 3) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 3) * a_chief * 1e3)
hold off
grid on
ylabel("\deltae_x' [m]")

subplot(3,2,2)
hold on
plot(t / T, mod_ROE_gt(:, 4) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 4) * a_chief * 1e3)
hold off
grid on
ylabel("\deltae_y' [m]")

subplot(3,2,4)
hold on
plot(t / T, mod_ROE_gt(:, 5) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 5) * a_chief * 1e3)
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(t / T, mod_ROE_gt(:, 6) * a_chief * 1e3)
plot(t / T, mod_ROE_stm(:, 6) * a_chief * 1e3)
hold off
grid on
ylabel('\deltai_y [m]')

figure
subplot(3,1,1)
hold on
plot(mod_ROE_gt(:, 2), mod_ROE_gt(:, 1))
plot(mod_ROE_stm(:, 2), mod_ROE_stm(:, 1))
hold off
axis equal
grid on
xlabel('\delta\lambda [m]')
ylabel('\deltaa [m]')

subplot(3,1,2)
hold on
plot(mod_ROE_gt(:,3), mod_ROE_gt(:,4))
plot(mod_ROE_stm(:,3), mod_ROE_stm(:,4))
hold off
axis equal
grid on
xlabel("\deltae_x' [m]")
ylabel('\deltae_y" [m]')

subplot(3,1,3)
hold on
plot(mod_ROE_gt(:,5), mod_ROE_gt(:,6))
plot(mod_ROE_stm(:,5), mod_ROE_stm(:,6))
hold off
axis equal
grid on
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')

figure
subplot(2,2,1)
plot(rel_motion_RTN(:,2), rel_motion_RTN(:,1))
axis equal
grid on
xlabel('T-axis [km]')
ylabel('R-axis [km]')

subplot(2,2,2)
plot(rel_motion_RTN(:,3), rel_motion_RTN(:,1))
axis equal
grid on
xlabel('N-axis [km]')
ylabel('R-axis [km]')

subplot(2,2,3)
plot(rel_motion_RTN(:,2), rel_motion_RTN(:,3))
axis equal
grid on
xlabel('T-axis [km]')
ylabel('N-axis [km]')

subplot(2,2,4)
plot3(rel_motion_RTN(:,1), rel_motion_RTN(:,2), rel_motion_RTN(:,3))
view(3)
axis equal
grid on
xlabel('R-axis [km]')
ylabel('T-axis [km]')
zlabel('N-axis [km]')

figure
subplot(2,2,1)
plot(rel_motion_RTN(:,5), rel_motion_RTN(:,4))
axis equal
grid on
xlabel('T-axis [km/s]')
ylabel('R-axis [km/s]')

subplot(2,2,2)
plot(rel_motion_RTN(:,6), rel_motion_RTN(:,4))
axis equal
grid on
xlabel('N-axis [km/s]')
ylabel('R-axis [km/s]')

subplot(2,2,3)
plot(rel_motion_RTN(:,5), rel_motion_RTN(:,6))
axis equal
grid on
xlabel('T-axis [km/s]')
ylabel('N-axis [km/s]')

subplot(2,2,4)
plot3(rel_motion_RTN(:,4), rel_motion_RTN(:,5), rel_motion_RTN(:,6))
view(3)
axis equal
grid on
xlabel('R-axis [km/s]')
ylabel('T-axis [km/s]')
zlabel('N-axis [km/S]')


%% Reconfiguration from Perigee Pass Mode to Inertial Attitude Mode using least squares

% Defining the ROE for Inertial Attitude Mode (IAM)

% Initial conditions for the chief
M_chief = n_chief * (6 * 3600 + 49 * 60);
nu_chief = mean2true(M_chief, e_chief);
oe_chief_IM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

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
nu_deputy = mean2true(M_deputy, e_deputy);
oe_deputy_IAM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

init_mod_ROE_IAM = OE2ROE_modified(oe_chief_PPM, oe_deputy_PPM)';



% % Initial conditions for the chief at the end of the perigee pass
% 
% a_chief = 36943; % km
% e_chief = 0.0001;
% inc_chief = deg2rad(59);
% omega_chief = deg2rad(188);
% RAAN_chief = deg2rad(84);
% n_chief = sqrt(mu / a_chief^3);
% T = 2*pi/n_chief;
% M_chief = n_chief * (4 * 3600 + 49 * 60);
% nu_chief = mean2true(M_chief, e_chief);
% % oe_chief = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];
% 
% % ROEs at the end of the perigee pass
% delta_a = 0;
% delta_lambda = 0 / a_chief;
% delta_ex = 800e-3 / a_chief;
% delta_ey = 200e-3 / a_chief;
% delta_ix = 750e-3 / a_chief;
% delta_iy = 200e-3 / a_chief;
% init_ROE = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]';
% 
% a_deputy = a_chief + delta_a * a_chief;
% inc_deputy = inc_chief + delta_ix;
% RAAN_deputy = RAAN_chief + delta_iy / sin(inc_chief);
% e_deputy = sqrt((delta_ex + e_chief * cos(omega_chief))^2 + (delta_ey + e_chief * sin(omega_chief))^2);
% omega_deputy = atan2(delta_ey + e_chief * sin(omega_chief), delta_ex + e_chief * cos(omega_chief));
% M_deputy = delta_lambda + M_chief + omega_chief - (RAAN_deputy - RAAN_chief) * cos(inc_chief) - omega_deputy;
% nu_deputy = mean2true(M_deputy, e_deputy);
% % oe_deputy = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];
% 
% % Position and velocity in ECI
% [pos_chief, vel_chief] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief, mu);
% [pos_deputy, vel_deputy] = OE2ECI(a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, nu_deputy, mu);
% 
% final_ROE = [0, 0, 59e-3/a_chief, 30e-3/a_chief, 80e-3/a_chief, 30e-3/a_chief]';
% STM_f0 = STM_ROE(2 * 3600, n_chief);
% delta_alpha = final_ROE - STM_f0 * init_ROE;
% 
% % Computing the maneuvers for the least squares
% STM_0 = eye(6);
% STM_1 = STM_ROE(30 * 60, n_chief);
% STM_2 = STM_ROE(3600, n_chief);
% STM_3 = STM_ROE(1.5 * 3600, n_chief);
% 
% B_0 = control_matrix(nu_chief, e_chief, inc_chief, a_chief, n_chief, omega_chief);
% B_1 = control_matrix(mean2true(M_chief + 30 * 60 * n_chief, e_chief), e_chief, inc_chief, a_chief, n_chief, omega_chief);
% B_2 = control_matrix(mean2true(M_chief + 3600 * n_chief, e_chief), e_chief, inc_chief, a_chief, n_chief, omega_chief);
% B_3 = control_matrix(mean2true(M_chief + 1.5 * 3600 * n_chief, e_chief), e_chief, inc_chief, a_chief, n_chief, omega_chief);
% B_4 = control_matrix(mean2true(M_chief + 2 * 3600 * n_chief, e_chief), e_chief, inc_chief, a_chief, n_chief, omega_chief);
% 
% super_control = [STM_f0 * B_0, STM_3 * B_1, STM_2 * B_2, STM_1 * B_3, STM_0 * B_4];
% delta_v = (super_control' * super_control)^(-1) * super_control' * delta_alpha; % in km/s
% 
% % Apply the first maneuver to the initial conditions
% 
% state_deputy_RTN = ECI2RTN([pos_chief; vel_chief], [pos_deputy; vel_deputy], mu);
% state_deputy_RTN(4:6) = state_deputy_RTN(4:6) + delta_v(1:3)';
% state_deputy_ECI = RTN2ECI([pos_chief; vel_chief], state_deputy_RTN, mu);
% 
% tspan_01 = linspace(4 * 3600 + 49 * 60, 4 * 3600 + 49 * 60 + 30 * 60, 1000); % simulating for half an hour between the first and second maneuver
% init_cond_01 = [pos_chief', vel_chief', state_deputy_ECI];
% 
% % Numerical integration without J2 perturbations
% [t_01, y_01] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_01, init_cond_01, options);
% 
% ROE_gt_01 = zeros(1000,6);
% 
% for j=1:1000
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_01(j,1:6), mu);
%     oe_chief = [a, e, i, omega, RAAN, M];
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_01(j,7:12), mu);
%     oe_deputy = [a, e, i, omega, RAAN, M];
%     ROE = OE2ROE(oe_chief, oe_deputy);
%     ROE_gt_01(j,:) = ROE;
% end
% 
% state_deputy_RTN_1 = ECI2RTN(y_01(end, 1:6), y_01(end, 7:12), mu);
% state_deputy_RTN_1(4:6) = state_deputy_RTN_1(4:6) + delta_v(4:6)';
% state_deputy_ECI_1 = RTN2ECI(y_01(end, 1:6), state_deputy_RTN_1, mu);
% 
% init_cond_12 = [y_01(end, 1:6), state_deputy_ECI_1];
% tspan_12 = linspace(4 * 3600 + 49 * 60 + 30 * 60, 5 * 3600 + 49 * 60, 1000);
% 
% % Numerical integration without J2 perturbations
% [t_12, y_12] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_12, init_cond_12, options);
% 
% ROE_gt_12 = zeros(1000,6);
% 
% for j=1:1000
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_12(j,1:6), mu);
%     oe_chief = [a, e, i, omega, RAAN, M];
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_12(j,7:12), mu);
%     oe_deputy = [a, e, i, omega, RAAN, M];
%     ROE = OE2ROE(oe_chief, oe_deputy);
%     ROE_gt_12(j,:) = ROE;
% end
% 
% state_deputy_RTN_2 = ECI2RTN(y_12(end, 1:6), y_12(end, 7:12), mu);
% state_deputy_RTN_2(4:6) = state_deputy_RTN_2(4:6) + delta_v(7:9)';
% state_deputy_ECI_2 = RTN2ECI(y_12(end, 1:6), state_deputy_RTN_2, mu);
% 
% init_cond_23 = [y_12(end, 1:6), state_deputy_ECI_2];
% tspan_23 = linspace(5 * 3600 + 49 * 60, 5 * 3600 + 79 * 60, 1000);
% 
% % Numerical integration without J2 perturbations
% [t_23, y_23] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_23, init_cond_23, options);
% 
% ROE_gt_23 = zeros(1000,6);
% 
% for j=1:1000
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_23(j,1:6), mu);
%     oe_chief = [a, e, i, omega, RAAN, M];
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_23(j,7:12), mu);
%     oe_deputy = [a, e, i, omega, RAAN, M];
%     ROE = OE2ROE(oe_chief, oe_deputy);
%     ROE_gt_23(j,:) = ROE;
% end
% 
% state_deputy_RTN_3 = ECI2RTN(y_23(end, 1:6), y_23(end, 7:12), mu);
% state_deputy_RTN_3(4:6) = state_deputy_RTN_3(4:6) + delta_v(10:12)';
% state_deputy_ECI_3 = RTN2ECI(y_23(end, 1:6), state_deputy_RTN_3, mu);
% 
% init_cond_34 = [y_23(end, 1:6), state_deputy_ECI_3];
% tspan_34 = linspace(5 * 3600 + 79 * 60, 6 * 3600 + 49 * 60, 1000);
% 
% % Numerical integration without J2 perturbations
% [t_34, y_34] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_34, init_cond_34, options);
% 
% ROE_gt_34 = zeros(1000,6);
% 
% for j=1:1000
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_34(j,1:6), mu);
%     oe_chief = [a, e, i, omega, RAAN, M];
%     [a, e, i, omega, RAAN, M] = ECI2OE_M(y_34(j,7:12), mu);
%     oe_deputy = [a, e, i, omega, RAAN, M];
%     ROE = OE2ROE(oe_chief, oe_deputy);
%     ROE_gt_34(j,:) = ROE;
% end
% 
% state_deputy_RTN_4 = ECI2RTN(y_34(end, 1:6), y_34(end, 7:12), mu);
% state_deputy_RTN_4(4:6) = state_deputy_RTN_4(4:6) + delta_v(13:15)';
% state_deputy_ECI_4 = RTN2ECI(y_34(end, 1:6), state_deputy_RTN_4, mu);
% 
% [a, e, i, omega, RAAN, M] = ECI2OE_M(y_34(end,1:6), mu);
% oe_chief = [a, e, i, omega, RAAN, M];
% [a, e, i, omega, RAAN, M] = ECI2OE_M(state_deputy_ECI_4, mu);
% oe_deputy = [a, e, i, omega, RAAN, M];
% ROE_final = OE2ROE(oe_chief, oe_deputy);
% 
% figure
% subplot(3,2,1)
% hold on
% plot((4 * 3600 + 49 * 60) / T, delta_a * a_chief * 1e3, 'o')
% plot(t_01 / T, ROE_gt_01(:, 1) * a_chief * 1e3)
% plot(t_12 / T, ROE_gt_12(:, 1) * a_chief * 1e3)
% plot(t_23 / T, ROE_gt_23(:, 1) * a_chief * 1e3)
% plot(t_34 / T, ROE_gt_34(:, 1) * a_chief * 1e3)
% plot(t_34(end) / T, ROE_final(1) * a_chief * 1e3, 'x')
% hold off
% grid on
% % legend('Ground Truth', 'STM')
% ylabel('\deltaa')
% 
% subplot(3,2,3)
% hold on
% plot((4 * 3600 + 49 * 60) / T, delta_lambda * a_chief * 1e3, 'o')
% plot(t_01 / T, ROE_gt_01(:, 2) * a_chief * 1e3)
% plot(t_12 / T, ROE_gt_12(:, 2) * a_chief * 1e3)
% plot(t_23 / T, ROE_gt_23(:, 2) * a_chief * 1e3)
% plot(t_34 / T, ROE_gt_34(:, 2) * a_chief * 1e3)
% plot(t_34(end) / T, ROE_final(2) * a_chief * 1e3, 'x')
% hold off
% grid on
% ylabel('\delta\lambda')
% 
% subplot(3,2,5)
% hold on
% plot((4 * 3600 + 49 * 60) / T, delta_ex * a_chief * 1e3, 'o')
% plot(t_01 / T, ROE_gt_01(:, 3) * a_chief * 1e3)
% plot(t_12 / T, ROE_gt_12(:, 3) * a_chief * 1e3)
% plot(t_23 / T, ROE_gt_23(:, 3) * a_chief * 1e3)
% plot(t_34 / T, ROE_gt_34(:, 3) * a_chief * 1e3)
% plot(t_34(end) / T, ROE_final(3) * a_chief * 1e3, 'x')
% hold off
% grid on
% ylabel('\deltae_x')
% 
% subplot(3,2,2)
% hold on
% plot((4 * 3600 + 49 * 60) / T, delta_ey * a_chief * 1e3, 'o')
% plot(t_01 / T, ROE_gt_01(:, 4) * a_chief * 1e3)
% plot(t_12 / T, ROE_gt_12(:, 4) * a_chief * 1e3)
% plot(t_23 / T, ROE_gt_23(:, 4) * a_chief * 1e3)
% plot(t_34 / T, ROE_gt_34(:, 4) * a_chief * 1e3)
% plot(t_34(end) / T, ROE_final(4) * a_chief * 1e3, 'x')
% hold off
% grid on
% ylabel('\deltae_y')
% 
% subplot(3,2,4)
% hold on
% plot((4 * 3600 + 49 * 60) / T, delta_ix * a_chief * 1e3, 'o')
% plot(t_01 / T, ROE_gt_01(:, 5) * a_chief * 1e3)
% plot(t_12 / T, ROE_gt_12(:, 5) * a_chief * 1e3)
% plot(t_23 / T, ROE_gt_23(:, 5) * a_chief * 1e3)
% plot(t_34 / T, ROE_gt_34(:, 5) * a_chief * 1e3)
% plot(t_34(end) / T, ROE_final(5) * a_chief * 1e3, 'x')
% hold off
% grid on
% ylabel('\deltai_x')
% 
% subplot(3,2,6)
% hold on
% plot((4 * 3600 + 49 * 60) / T, delta_iy * a_chief * 1e3, 'o')
% plot(t_01 / T, ROE_gt_01(:, 6) * a_chief * 1e3)
% plot(t_12 / T, ROE_gt_12(:, 6) * a_chief * 1e3)
% plot(t_23 / T, ROE_gt_23(:, 6) * a_chief * 1e3)
% plot(t_34 / T, ROE_gt_34(:, 6) * a_chief * 1e3)
% plot(t_34(end) / T, ROE_final(6) * a_chief * 1e3, 'x')
% hold off
% grid on
% ylabel('\deltai_y')
% 
% % Then simulate until the time of each maneuver, compute the new initial
% % conditions and keep propagating
%                                                                                                                   


%% Funtions

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end

function ROE = OE2ROE_modified(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) - oe_chief(2);
    ROE(4) = wrapToPi(oe_deputy(4) - oe_chief(4) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end

function STM = STM_ROE(t, n)
    STM = eye(6);
    STM(2,1) = -1.5 * n * t;
end

function B = control_matrix(nu, e, i, a, n, omega)
    eta = sqrt(1-e^2);                   
    B = zeros(6, 3);
    B(1,1) = 2 * e * sin(nu) / eta;
    B(1,2) = (1+e*cos(nu)) * 2 / eta;
    B(2,1) = - 2 * eta^2 / (1+e*cos(nu));
    B(3,1) = eta * sin(nu);
    B(3,2) = ((2+e*cos(nu))*cos(nu) + e) / (1+e*cos(nu)) * eta;
%     B(3,3) = eta * e * sin(omega) * sin(nu+omega) / ((tan(i)) * (1+e*cos(nu)));
    B(4,1) = - eta * cos(nu) / e;
    B(4,2) = ((2+e*cos(nu))*sin(nu)) / (1+e*cos(nu)) * eta / e;
%     B(4,3) = - eta * e * cos(omega) * sin(nu+omega) / ((tan(i)) * (1+e*cos(nu)));
    B(5,3) = eta * cos(nu+omega) / (1+e*cos(nu));
    B(6,3) = eta * sin(nu+omega) / (1+e*cos(nu));
    B = B / (n * a);
end

function B = control_matrix_2(nu, e, i, a, n, omega)
    eta = sqrt(1-e^2);                   
    B = zeros(6, 3);
    B(1,1) = 2 * e * sin(nu) / eta;
    B(1,2) = (1+e*cos(nu)) * 2 / eta;
    B(2,1) = - 2 * eta^2 / (1+e*cos(nu));
    B(3,1) = eta * sin(nu + omega);
    B(3,2) = ((2+e*cos(nu))*cos(nu+omega) + e*cos(omega)) / (1+e*cos(nu)) * eta;
    B(3,3) = eta * e * sin(omega) * sin(nu+omega) / ((tan(i)) * (1+e*cos(nu)));
    B(4,1) = - eta * cos(nu+omega);
    B(4,2) = ((2+e*cos(nu))*sin(nu+omega) + e*sin(omega)) / (1+e*cos(nu)) * eta;
    B(4,3) = - eta * e * cos(omega) * sin(nu+omega) / ((tan(i)) * (1+e*cos(nu)));
    B(5,3) = eta * cos(nu+omega) / (1+e*cos(nu));
    B(6,3) = eta * sin(nu+omega) / (1+e*cos(nu));
    B = B / (n * a);
end

function M = STM_ROE_2(f, e, omega, i, t, n)
    k = 1 + e * cos(f);
    k_prime = - e * sin(f);
    eta = sqrt(1 - e^2);
    u = omega + f;
    e_x = e * cos(omega);
    e_y = e * sin(omega);

    M = zeros(6);

    M(1, 1) = 1/k + 3/2 * k_prime * n / eta^3 * t;
    M(1, 2) = - k_prime / eta^3;
    M(1, 3) = 1/eta^3 * (e_x * (k-1)/(1+eta) - cos(u));
    M(1, 4) = 1/eta^3 * (e_y * (k-1)/(1+eta) - sin(u));
    M(1, 6) = k_prime / eta^3 * cot(i);

    M(2, 1) = -3/2 * k * n / eta^3 * t;
    M(2, 2) = k / eta^3;
    M(2, 3) = 1/eta^2 * ((1+1/k) * sin(u) + e_y/k + k/eta * e_y /(1+eta));
    M(2, 4) = - 1/eta^2 * ((1+1/k) * cos(u) + e_x/k + k/eta * e_x /(1+eta));
    M(2, 6) = (1/k - k/eta^3) * cot(i);

    M(3, 5) = 1/k * sin(u);
    M(3, 6) = -1/k * cos(u);

    M(4, 1) = k_prime/2 + 3/2 * k^2 * (1-k) * n /eta^3 * t;
    M(4, 2) = k^2 / eta^3 * (k-1);
    M(4, 3) = k^2/eta^3 * (eta*sin(u) + e_y * (k-1)/(1+eta));
    M(4, 4) = - k^2/eta^3 * (eta*cos(u) + e_x * (k-1)/(1+eta));
    M(4, 6) = -k^2/eta^3 * (k-1) * cot(i);

    M(5, 1) = -3/2 * k * (1 + k * k_prime * n / eta^3 * t);
    M(5, 2) = k^2 / eta^3 * k_prime;
    M(5, 3) = (1+k^2/eta^3)*cos(u) + e_x * k/eta^2 * (1 + k/eta * (1-k)/(1+eta));
    M(5, 4) = (1+k^2/eta^3)*sin(u) + e_y * k/eta^2 * (1 + k/eta * (1-k)/(1+eta));
    M(5, 6) = -(1+k^2/eta^3) * k_prime * cot(i);

    M(6, 5) = cos(u) + e_x;
    M(6, 6) = sin(u) + e_y;
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