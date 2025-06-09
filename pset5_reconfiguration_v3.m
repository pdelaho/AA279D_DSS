%% Reconfiguration problem

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100); % options for numerical integration

% Initial configuration (Perigee Pass Mode PPM)

a_chief = 36943; % km
e_chief = 0.8111;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
T_chief = 2 * pi / n_chief;
M_chief = n_chief * (4 * 3600 + 49 * 60);
nu_chief_PPM = mean2true(M_chief, e_chief, 1e-12);
oe_chief_PPM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

delta_a = 0;
delta_lambda = 100e-3 / a_chief;
delta_ex = 750e-3 / a_chief;
delta_ey = 150e-3 / a_chief;
delta_ix = 700e-3 / a_chief;
delta_iy = 140e-3 / a_chief;
% delta_a = 0e-3 / a_chief;
% delta_lambda = 0e-3 / a_chief;
% delta_ex = 0e-3 / a_chief;
% delta_ey = 0e-3 / a_chief;
% delta_ix = 0e-3 / a_chief;
% delta_iy = 0e-3 / a_chief;
ROE_PPM = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];

oe_deputy_PPM = ROE2OE(oe_chief_PPM, ROE_PPM);
ecc_ROE_PPM = eccentric_OE2ROE(oe_chief_PPM, oe_deputy_PPM);

% Final configuration mode (Inertial Attitude Mode IAM)

a_chief = 36943; % km
e_chief = 0.8111;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
% T_chief = 2 * pi / n_chief;
M_chief = n_chief * (6 * 3600 + 49 * 60);
% nu_chief_IAM = mean2true(M_chief, e_chief, 1e-12);
oe_chief_IAM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

delta_a = 0;
delta_lambda = 80e-3 / a_chief;
delta_ex = -20e-3 / a_chief;
delta_ey = 89.4e-3 / a_chief;
delta_ix = 0 / a_chief;
delta_iy = 79e-3 / a_chief;
% delta_a = 100e-3 / a_chief;
% delta_lambda = 0e-3 / a_chief;
% delta_ex = 0e-3 / a_chief;
% delta_ey = 0e-3 / a_chief;
% delta_ix = 0e-3 / a_chief;
% delta_iy = 0e-3 / a_chief;
ROE_IAM = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];

oe_deputy_IAM = ROE2OE(oe_chief_IAM, ROE_IAM);
ecc_ROE_IAM = eccentric_OE2ROE(oe_chief_IAM, oe_deputy_IAM);

% Computing impulsive maneuvers
% A_0f = STM(n_chief, 2 * 60 * 60);
% delta_ROE = ecc_ROE_IAM - A_0f * ecc_ROE_PPM;
% 
% B_0 = control_input(a_chief, n_chief, e_chief, nu_chief_PPM, omega_chief, inc_chief);
% 
% A_1f = STM(n_chief, 3600 + 40 * 60);
% B_1 = control_input(a_chief, n_chief, e_chief, mean2true((5 * 3600 + 9 * 60) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);
% 
% A_2f = STM(n_chief, 3600 + 20 * 60);
% B_2 = control_input(a_chief, n_chief, e_chief, mean2true((5 * 3600 + 29 * 60) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);
% 
% A_3f = STM(n_chief, 3600);
% B_3 = control_input(a_chief, n_chief, e_chief, mean2true((5 * 3600 + 49 * 60) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);
% 
% A_4f = STM(n_chief, 40 * 60);
% B_4 = control_input(a_chief, n_chief, e_chief, mean2true((6 * 3600 + 9 * 60) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);
% 
% A_5f = STM(n_chief, 20 * 60);
% B_5 = control_input(a_chief, n_chief, e_chief, mean2true((6 * 3600 + 29 * 60) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);

maneuvers_T = [T_chief, 5 * T_chief / 6, 4 * T_chief / 6, 3 * T_chief / 6, 2 * T_chief / 6, T_chief / 6];

A_0f = STM(n_chief, T_chief);
delta_ROE = ecc_ROE_IAM - A_0f * ecc_ROE_PPM;

B_0 = control_input(a_chief, n_chief, e_chief, nu_chief_PPM, omega_chief, inc_chief);

A_1f = STM(n_chief, 5 * T_chief / 6);
B_1 = control_input(a_chief, n_chief, e_chief, mean2true(5 * T_chief / 6, e_chief, 1e-14), omega_chief, inc_chief);

A_2f = STM(n_chief, 4 * T_chief / 6);
B_2 = control_input(a_chief, n_chief, e_chief, mean2true(4 * T_chief / 6, e_chief, 1e-14), omega_chief, inc_chief);

A_3f = STM(n_chief, 3 * T_chief / 6);
B_3 = control_input(a_chief, n_chief, e_chief, mean2true(3 * T_chief / 6, e_chief, 1e-14), omega_chief, inc_chief);

A_4f = STM(n_chief, 2 * T_chief / 6);
B_4 = control_input(a_chief, n_chief, e_chief, mean2true(2 * T_chief / 6, e_chief, 1e-14), omega_chief, inc_chief);

A_5f = STM(n_chief, T_chief / 6);
B_5 = control_input(a_chief, n_chief, e_chief, mean2true(T_chief / 6, e_chief, 1e-14), omega_chief, inc_chief);

G = [A_0f * B_0, A_1f * B_1, A_2f * B_2, A_3f * B_3, A_4f * B_4, A_5f * B_5];
delta_v = pinv(G) * delta_ROE;

% Apply impulsive maneuvers
[pos_chief_0, vel_chief_0] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief_PPM, mu);
[pos_deputy_0, vel_deputy_0] = OE2ECI(oe_deputy_PPM(1), oe_deputy_PPM(2), oe_deputy_PPM(3), oe_deputy_PPM(4), oe_deputy_PPM(5), mean2true(oe_deputy_PPM(6), oe_deputy_PPM(2), 1e-14), mu);

R_vec = pos_chief_0 / norm(pos_chief_0);
h_vec = cross(pos_chief_0, vel_chief_0);
N_vec = h_vec / norm(h_vec);
T_vec = cross(N_vec, R_vec);
rotation = [R_vec'; T_vec'; N_vec'];

% R_vec = pos_deputy_0 / norm(pos_deputy_0);
% h_vec = cross(pos_deputy_0, vel_deputy_0);
% N_vec = h_vec / norm(h_vec);
% T_vec = cross(N_vec, R_vec);
% rotation = [R_vec'; T_vec'; N_vec'];

vel_deputy_0 = vel_deputy_0 + rotation' * delta_v(1:3);
init_cond_0 = [pos_chief_0', vel_chief_0', pos_deputy_0', vel_deputy_0'];
% tspan_01 = linspace(4 * 3600 + 49 * 60, 5 * 3600 + 9 * 60, 1000);
tspan_01 = linspace(0, maneuvers_T(end), 1000);
[t_01, y_01] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_01 , init_cond_0, options);

ROE_01 = zeros(1000, 6);
for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_01(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_01(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = OE2ROE(oe_chief, oe_deputy);
    ROE_01(j,:) = ROE;
end

% R_vec = y_01(end, 1:3) / norm(y_01(end, 1:3));
% h_vec = cross(y_01(end, 1:3), y_01(end, 4:6));
% N_vec = h_vec / norm(h_vec);
% T_vec = cross(N_vec, R_vec);
% rotation = [R_vec; T_vec; N_vec];

R_vec = y_01(end, 7:9) / norm(y_01(end, 7:9));
h_vec = cross(y_01(end, 7:9), y_01(end, 10:12));
N_vec = h_vec / norm(h_vec);
T_vec = cross(N_vec, R_vec);
rotation = [R_vec; T_vec; N_vec];

vel_deputy_1 = y_01(end, 10:12)' + rotation' * delta_v(4:6);
init_cond_1 = [y_01(end, 1:9), vel_deputy_1'];
% tspan_12 = linspace(5 * 3600 + 9 * 60, 5 * 3600 + 29 * 60, 1000);
tspan_12 = linspace(maneuvers_T(end), maneuvers_T(end-1), 1000);
[t_12, y_12] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_12 , init_cond_1, options);

ROE_12 = zeros(1000, 6);
for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_12(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_12(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = OE2ROE(oe_chief, oe_deputy);
    ROE_12(j,:) = ROE;
end

% R_vec = y_12(end, 1:3) / norm(y_12(end, 1:3));
% h_vec = cross(y_12(end, 1:3), y_12(end, 4:6));
% N_vec = h_vec / norm(h_vec);
% T_vec = cross(N_vec, R_vec);
% rotation = [R_vec; T_vec; N_vec];

R_vec = y_12(end, 7:9) / norm(y_12(end, 7:9));
h_vec = cross(y_12(end, 7:9), y_12(end, 10:12));
N_vec = h_vec / norm(h_vec);
T_vec = cross(N_vec, R_vec);
rotation = [R_vec; T_vec; N_vec];

vel_deputy_2 = y_12(end, 10:12)' + rotation' * delta_v(7:9);
init_cond_2 = [y_12(end, 1:9), vel_deputy_2'];
% tspan_23 = linspace(5 * 3600 + 29 * 60, 5 * 3600 + 49 * 60, 1000);
tspan_23 = linspace(maneuvers_T(end-1), maneuvers_T(end-2), 1000);
[t_23, y_23] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_23, init_cond_2, options);

ROE_23 = zeros(1000, 6);
for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_23(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_23(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = OE2ROE(oe_chief, oe_deputy);
    ROE_23(j,:) = ROE;
end

% R_vec = y_23(end, 1:3) / norm(y_23(end, 1:3));
% h_vec = cross(y_23(end, 1:3), y_23(end, 4:6));
% N_vec = h_vec / norm(h_vec);
% T_vec = cross(N_vec, R_vec);
% rotation = [R_vec; T_vec; N_vec];

R_vec = y_23(end, 7:9) / norm(y_23(end, 7:9));
h_vec = cross(y_23(end, 7:9), y_23(end, 10:12));
N_vec = h_vec / norm(h_vec);
T_vec = cross(N_vec, R_vec);
rotation = [R_vec; T_vec; N_vec];

vel_deputy_3 = y_23(end, 10:12)' + rotation' * delta_v(10:12);
init_cond_3 = [y_23(end, 1:9), vel_deputy_3'];
% tspan_34 = linspace(5 * 3600 + 49 * 60, 6 * 3600 + 9 * 60, 1000);
tspan_34 = linspace(maneuvers_T(end-2), maneuvers_T(end-3), 1000);
[t_34, y_34] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_34 , init_cond_3, options);

ROE_34 = zeros(1000, 6);
for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_34(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_34(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = OE2ROE(oe_chief, oe_deputy);
    ROE_34(j,:) = ROE;
end

R_vec = y_34(end, 1:3) / norm(y_34(end, 1:3));
h_vec = cross(y_34(end, 1:3), y_34(end, 4:6));
N_vec = h_vec / norm(h_vec);
T_vec = cross(N_vec, R_vec);
rotation = [R_vec; T_vec; N_vec];

% R_vec = y_34(end, 7:9) / norm(y_34(end, 7:9));
% h_vec = cross(y_34(end, 7:9), y_34(end, 7:9));
% N_vec = h_vec / norm(h_vec);
% T_vec = cross(N_vec, R_vec);
% rotation = [R_vec; T_vec; N_vec];

vel_deputy_4 = y_34(end, 10:12)' + rotation' * delta_v(13:15);
init_cond_4 = [y_34(end, 1:9), vel_deputy_4'];
% tspan_45 = linspace(6 * 3600 + 9 * 60, 6 * 3600 + 29 * 60, 1000);
tspan_45 = linspace(maneuvers_T(end-3), maneuvers_T(end-4), 1000);
[t_45, y_45] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_45 , init_cond_4, options);

ROE_45 = zeros(1000, 6);
for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_45(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_45(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = OE2ROE(oe_chief, oe_deputy);
    ROE_45(j,:) = ROE;
end

% R_vec = y_45(end, 1:3) / norm(y_45(end, 1:3));
% h_vec = cross(y_45(end, 1:3), y_45(end, 4:6));
% N_vec = h_vec / norm(h_vec);
% T_vec = cross(N_vec, R_vec);
% rotation = [R_vec; T_vec; N_vec];

R_vec = y_45(end, 7:9) / norm(y_45(end, 7:9));
h_vec = cross(y_45(end, 7:9), y_45(end, 10:12));
N_vec = h_vec / norm(h_vec);
T_vec = cross(N_vec, R_vec);
rotation = [R_vec; T_vec; N_vec];

vel_deputy_5 = y_45(end, 10:12)' + rotation' * delta_v(16:18);
init_cond_5 = [y_45(end, 1:9), vel_deputy_5'];
% tspan_56 = linspace(6 * 3600 + 29 * 60, 6 * 3600 + 49 * 60, 1000);
tspan_56 = linspace(maneuvers_T(end-4), maneuvers_T(end-5), 1000);
[t_56, y_56] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan_56, init_cond_5, options);

ROE_56 = zeros(1000, 6);
for j=1:1000
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_56(j,1:6), mu);
    oe_chief = [a, e, i, omega, RAAN, M];
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_56(j,7:12), mu);
    oe_deputy = [a, e, i, omega, RAAN, M];
    ROE = OE2ROE(oe_chief, oe_deputy);
    ROE_56(j,:) = ROE;
end

%% Figures

figure
subplot(3,2,1)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, ROE_PPM(1) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, ROE_IAM(1) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, ROE_01(:, 1) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:, 1) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:, 1) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:, 1) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:, 1) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:, 1) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
legend('Initial', 'Target')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, ROE_PPM(2) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_01(:, 2) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:, 2) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:, 2) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:, 2) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:, 2) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:, 2) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, ROE_IAM(2) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, ROE_PPM(3) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_01(:, 3) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:, 3) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:, 3) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:, 3) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:, 3) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:, 3) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, ROE_IAM(3) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel("\deltae_x [m]")
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, ROE_PPM(4) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_01(:, 4) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:, 4) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:, 4) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:, 4) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:, 4) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:, 4) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, ROE_IAM(4) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
ylabel("\deltae_y [m]")

subplot(3,2,4)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, ROE_PPM(5) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_01(:, 5) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:, 5) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:, 5) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:, 5) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:, 5) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:, 5) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, ROE_IAM(5) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot((4 * 3600 + 49 * 60) / T_chief, ROE_PPM(6) * a_chief * 1e3, 'ro')
plot(t_01 / T_chief, ROE_01(:, 6) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:, 6) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:, 6) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:, 6) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:, 6) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:, 6) * a_chief * 1e3, 'b-')
plot(t_56(end) / T_chief, ROE_IAM(6) * a_chief * 1e3, 'rx')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltai_y [m]')
xlabel('Orbital Periods')

%% Functions

% function ROE = eccentric_OE2ROE(oe_chief, oe_deputy)
%     ROE = zeros(6, 1);
%     eta = sqrt(1 - oe_chief(2)^2);
% 
%     ROE(1) = eta^2 * (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
%     ROE(2) = (1 / eta) * (oe_deputy(6) - oe_chief(6)) + eta^2 * (oe_deputy(4) - oe_chief(4) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
%     ROE(3) = (oe_deputy(2) - oe_chief(2)) * cos(oe_chief(4)) + oe_chief(2) / eta * (oe_deputy(6) - oe_chief(6)) * sin(oe_chief(4));
%     ROE(4) = (oe_deputy(2) - oe_chief(2)) * sin(oe_chief(4)) - oe_chief(2) / eta * (oe_deputy(6) - oe_chief(6)) * cos(oe_chief(4));
%     ROE(5) = eta^2 * (oe_deputy(3) - oe_chief(3));
%     ROE(6) = eta^2 * (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
% end

% function ROE = eccentric_OE2ROE(oe_chief, oe_deputy)
%     ROE = zeros(6, 1);
%     eta = sqrt(1 - oe_chief(2)^2);
% 
%     ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
%     ROE(2) = (oe_deputy(6) - oe_chief(6)) + eta * (oe_deputy(4) - oe_chief(4) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
%     ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
%     ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
%     ROE(5) = oe_deputy(3) - oe_chief(3);
%     ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
% end

function ROE = eccentric_OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(6, 1);
    eta = sqrt(1 - oe_chief(2)^2);

    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = (oe_deputy(6) - oe_chief(6)) + eta * (oe_deputy(4) - oe_chief(4) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) - oe_chief(2);
    ROE(4) = oe_deputy(4) - oe_chief(4) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end

% function ROE = OE2ROE(oe_chief, oe_deputy)
%     ROE = zeros(1, 6);
%     ROE(1) = (oe_deputy(1) - oe_chief(1));
%     ROE(2) = wrapTo2Pi(((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3))) * oe_chief(1));
%     ROE(3) = (oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4))) * oe_chief(1);
%     ROE(4) = (oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4))) * oe_chief(1);
%     ROE(5) = (oe_deputy(3) - oe_chief(3)) * oe_chief(1);
%     ROE(6) = ((oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3))) * oe_chief(1);
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
    oe_deputy(1) = oe_chief(1) + ROE(1);
    oe_deputy(3) = oe_chief(3) + ROE(5) / oe_chief(1);
    oe_deputy(2) = sqrt((ROE(3) / oe_chief(1) + oe_chief(2) * cos(oe_chief(4)))^2 + (ROE(4) / oe_chief(1) + oe_chief(2) * sin(oe_chief(4)))^2);
    oe_deputy(5) = ROE(6) / oe_chief(1) / sin(oe_chief(3)) + oe_chief(5);
    oe_deputy(4) = atan2(ROE(4) / oe_chief(1) + oe_chief(2) * sin(oe_chief(4)), ROE(3) / oe_chief(1) + oe_chief(2) * cos(oe_chief(4)));
    oe_deputy(6) = ROE(2) / oe_chief(1) + oe_chief(6) - oe_deputy(4) + oe_chief(4) - (oe_deputy(5) - oe_chief(5))*cos(oe_chief(3));
end

% function B = control_input(a, n, e, nu, omega, i)
%     B = zeros(6, 3);
%     eta = sqrt(1 - e^2);
% 
%     B(1,1) = 2 * e * sin(nu);
%     B(1,2) = 2 * (1 + e *cos(nu));
% 
%     B(2,1) = (e * cos(nu) + 2) * (e * cos(nu) - 1) / (1 + e * cos(nu));
%     B(2,2) = (2 + e * cos(nu)) * (sin(nu) - eta^2) / (e * (1 + e * cos(nu)));
% 
%     B(3,1) = (1 + e * cos(nu) * sin(omega + nu) - 2 * e * sin(omega)) / (1 + e *cos(nu));
%     B(3,2) = ((2 + e * cos(nu)) * cos(omega + nu) + e * cos(omega)) / (1 + e * cos(nu));
% 
%     B(4,1) = -((1 + e * cos(nu)) * cos(omega + nu) - 2 * e * cos(omega)) / (1 + e * cos(nu));
%     B(4,2) = ((2 + e * cos(nu)) * sin(omega + nu) + e * sin(omega)) / (1 + e * cos(nu));
% 
%     B(5,3) = eta^2 * cos(omega + nu) / (1 + e *cos(nu));
% 
%     B(6,3) = eta^2 * sin(omega + nu) / (1 + e * cos(nu));
% 
%     B = B .* eta ./ n ./ a;
% end

% function B = control_input(a, n, e, nu, omega, i)
%     B = zeros(6, 3);
%     eta = sqrt(1 - e^2);
% 
%     B(1,1) = 2 * e * sin(nu) / eta;
%     B(1,2) = 2 * (1 + e *cos(nu)) / eta;
% 
%     B(2,1) = -2 * eta^2 / (1 + e * cos(nu));
% 
%     B(3,1) = eta * sin(omega + nu);
%     B(3,2) = ((2 + e * cos(nu)) * cos(omega + nu) + e * cos(omega)) / (1 + e * cos(nu)) * eta;
%     B(3,3) = eta * e * sin(omega) / tan(i) * sin(omega + nu) / (1 + e * cos(nu));
% 
%     B(4,1) = -eta * cos(omega + nu);
%     B(4,2) = ((2 + e * cos(nu)) * sin(omega + nu) + e * sin(omega)) / (1 + e * cos(nu)) * eta;
%     B(4,3) = -eta * e * cos(omega) / tan(i) * sin(omega + nu) / (1 + e * cos(nu));
% 
%     B(5,3) = eta * cos(omega + nu) / (1 + e *cos(nu));
% 
%     B(6,3) = eta * sin(omega + nu) / (1 + e * cos(nu));
% 
%     B = B ./ n ./ a;
% end

function B = control_input(a, n, e, nu, omega, i)
    B = zeros(6, 3);
    eta = sqrt(1 - e^2);

    B(1,1) = 2 * e * sin(nu) / eta;
    B(1,2) = 2 * (1 + e *cos(nu)) / eta;

    B(2,1) = -2 * eta^2 / (1 + e * cos(nu));

    B(3,1) = eta * sin(nu);
    B(3,2) = eta * (e + cos(nu) * (2 + e * cos(nu))) / (1 + e * cos(nu));

    B(4,1) = -eta * cos(nu) / e;
    B(4,2) = (eta / e) * sin(nu) * (2 + e * cos(nu)) / (1 + e * cos(nu));

    B(5,3) = eta * cos(omega + nu) / (1 + e *cos(nu));

    B(6,3) = eta * sin(omega + nu) / (1 + e * cos(nu));

    B = B ./ n ./ a;
end

function A = STM(n, delta_t)
    A = eye(6);
    A(2,1) = -1.5 * n * delta_t;
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