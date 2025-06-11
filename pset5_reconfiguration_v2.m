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

oe_deputy_PPM = ROE2OE(oe_chief_PPM, ROE_PPM)';
% ecc_ROE_PPM = eccentric_OE2ROE(oe_chief_PPM, oe_deputy_PPM);

% Final configuration mode (Inertial Attitude Mode IAM)

a_chief = 36943; % km
e_chief = 0.8111;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
% T_chief = 2 * pi / n_chief;
M_chief = n_chief * (6 * 3600 + 49 * 60);
nu_chief_IAM = mean2true(M_chief, e_chief, 1e-12);
oe_chief_IAM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

delta_a = 0;
delta_lambda = 80e-3 / a_chief;
delta_ex = -20e-3 / a_chief;
delta_ey = 89.4e-3 / a_chief;
delta_ix = 0 / a_chief;
delta_iy = 79e-3 / a_chief;
% delta_a = 0e-3 / a_chief;
% delta_lambda = 0e-3 / a_chief;
% delta_ex = 0e-3 / a_chief;
% delta_ey = 100000e-3 / a_chief;
% delta_ix = 0e-3 / a_chief;
% delta_iy = 0e-3 / a_chief;
ROE_IAM = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy];

oe_deputy_IAM = ROE2OE(oe_chief_IAM, ROE_IAM)';
% ecc_ROE_IAM = eccentric_OE2ROE(oe_chief_IAM, oe_deputy_IAM);

[pos_chief_IAM, vel_chief_IAM] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief_IAM, mu);
[pos_deputy_IAM, vel_deputy_IAM] = OE2ECI(oe_deputy_IAM(1), oe_deputy_IAM(2), oe_deputy_IAM(3), oe_deputy_IAM(4), oe_deputy_IAM(5), mean2true(oe_deputy_IAM(6), oe_deputy_IAM(2), 1e-12), mu);
rel_pos_IAM = pos_deputy_IAM - pos_chief_IAM;
rel_vel_IAM = vel_deputy_IAM - vel_chief_IAM;

% Computing delta-v
times = [4 * 3600 + 49 * 60, 5 * 3600 + 9 * 60, 5 * 3600 + 29 * 60, 5 * 3600 + 49 * 60, 6 * 3600 + 9 * 60, 6 * 3600 + 29 * 60, 6 * 3600 + 49 * 60];
A_0f = ROE_STM(n_chief, times(end) - times(1));
delta_ROE = ROE_IAM' - A_0f * ROE_PPM';

B_0 = control_input(a_chief, n_chief, e_chief, nu_chief_PPM, omega_chief, inc_chief);

A_1f = ROE_STM(n_chief, times(end) - times(2));
B_1 = control_input(a_chief, n_chief, e_chief, mean2true(times(2) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);

A_2f = ROE_STM(n_chief, times(end) - times(3));
B_2 = control_input(a_chief, n_chief, e_chief, mean2true(times(3) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);

A_3f = ROE_STM(n_chief, times(end) - times(4));
B_3 = control_input(a_chief, n_chief, e_chief, mean2true(times(4) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);

A_4f = ROE_STM(n_chief, times(end) - times(5));
B_4 = control_input(a_chief, n_chief, e_chief, mean2true(times(5) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);

A_5f = ROE_STM(n_chief, times(end) - times(6));
B_5 = control_input(a_chief, n_chief, e_chief, mean2true(times(6) * n_chief, e_chief, 1e-14), omega_chief, inc_chief);

G = [A_0f * B_0, A_1f * B_1, A_2f * B_2, A_3f * B_3, A_4f * B_4, A_5f * B_5];
delta_v = pinv(G) * delta_ROE;

% Applying delta-v
delta_ROE_0 = transform(oe_chief_PPM, delta_v(1:3), mu);
ic_01 = [oe_chief_PPM, oe_deputy_PPM + delta_ROE_0];
tspan_01 = linspace(times(1), times(2), 1000);
[t_01, y_01] = ode89(@ (t, state) GVE(t, state, mu), tspan_01, ic_01, options);


ROE_01 = zeros(1000, 6);
for j=1:1000
    ROE_01(j, :) = OE2ROE(y_01(j, 1:6), y_01(j, 7:12));
end

delta_ROE_1 = transform(y_01(end, 1:6), delta_v(4:6), mu);
ic_12 = [y_01(end, 1:6), y_01(end, 7:12) + delta_ROE_1];
tspan_12 = linspace(times(2), times(3), 1000);
[t_12, y_12] = ode89(@ (t, state) GVE(t, state, mu), tspan_12, ic_12, options);

ROE_12 = zeros(1000, 6);
for j=1:1000
    ROE_12(j, :) = OE2ROE(y_12(j, 1:6), y_12(j, 7:12));
end

delta_ROE_2 = transform(y_12(end, 1:6), delta_v(7:9), mu);
ic_23 = [y_12(end, 1:6), y_12(end, 7:12) + delta_ROE_2];
tspan_23 = linspace(times(3), times(4), 1000);
[t_23, y_23] = ode89(@ (t, state) GVE(t, state, mu), tspan_23, ic_23, options);

ROE_23 = zeros(1000, 6);
for j=1:1000
    ROE_23(j, :) = OE2ROE(y_23(j, 1:6), y_23(j, 7:12));
end

delta_ROE_3 = transform(y_23(end, 1:6), delta_v(10:12), mu);
ic_34 = [y_23(end, 1:6), y_23(end, 7:12) + delta_ROE_3];
tspan_34 = linspace(times(4), times(5), 1000);
[t_34, y_34] = ode89(@ (t, state) GVE(t, state, mu), tspan_34, ic_34, options);

ROE_34 = zeros(1000, 6);
for j=1:1000
    ROE_34(j, :) = OE2ROE(y_34(j, 1:6), y_34(j, 7:12));
end

delta_ROE_4 = transform(y_34(end, 1:6), delta_v(13:15), mu);
ic_45 = [y_34(end, 1:6), y_34(end, 7:12) + delta_ROE_4];
tspan_45 = linspace(times(5), times(6), 1000);
[t_45, y_45] = ode89(@ (t, state) GVE(t, state, mu), tspan_45, ic_45, options);

ROE_45 = zeros(1000, 6);
for j=1:1000
    ROE_45(j, :) = OE2ROE(y_45(j, 1:6), y_45(j, 7:12));
end

delta_ROE_5 = transform(y_45(end, 1:6), delta_v(16:18), mu);
ic_56 = [y_45(end, 1:6), y_45(end, 7:12) + delta_ROE_5];
tspan_56 = linspace(times(6), times(7), 1000);
[t_56, y_56] = ode89(@ (t, state) GVE(t, state, mu), tspan_56, ic_56, options);

ROE_56 = zeros(1000, 6);
for j=1:1000
    ROE_56(j, :) = OE2ROE(y_56(j, 1:6), y_56(j, 7:12));
end

[pos_chief_f, vel_chief_f] = OE2ECI(y_56(end, 1), y_56(end, 2), y_56(end, 3), y_56(end, 4), y_56(end, 5), mean2true(y_56(end, 6), y_56(end, 2), 1e-12), mu);
[pos_deputy_f, vel_deputy_f] = OE2ECI(y_56(end, 7), y_56(end, 8), y_56(end, 9), y_56(end, 10), y_56(end, 11), mean2true(y_56(end, 12), y_56(end, 8), 1e-12), mu);
rel_pos_f = pos_deputy_f - pos_chief_f;
rel_vel_f = vel_deputy_f - vel_chief_f;

%% Figures

figure
subplot(3,2,1)
hold on
plot(t_01(1) / T_chief, ROE_PPM(1) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, ROE_IAM(1) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, ROE_01(:,1) * a_chief * 1e3, 'b')
plot(t_12 / T_chief, ROE_12(:,1) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:,1) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:,1) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:,1) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:,1) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
legend('Start', 'Target', 'Maneuvering')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(t_01(1) / T_chief, ROE_PPM(2) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, ROE_IAM(2) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, ROE_01(:,2) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:,2) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:,2) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:,2) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:,2) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:,2) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(t_01(1) / T_chief, ROE_PPM(3) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, ROE_IAM(3) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, ROE_01(:,3) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:,3) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:,3) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:,3) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:,3) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:,3) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_01(1) / T_chief, ROE_PPM(4) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, ROE_IAM(4) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, ROE_01(:,4) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:,4) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:,4) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:,4) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:,4) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:,4) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(t_01(1) / T_chief, ROE_PPM(5) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, ROE_IAM(5) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, ROE_01(:,5) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:,5) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:,5) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:,5) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:,5) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:,5) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(t_01(1) / T_chief, ROE_PPM(6) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, ROE_IAM(6) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, ROE_01(:,6) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, ROE_12(:,6) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, ROE_23(:,6) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, ROE_34(:,6) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, ROE_45(:,6) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, ROE_56(:,6) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltaa [m]')
xlabel('Orbital Periods')


figure
subplot(3,2,1)
hold on
plot(t_01(1) / T_chief, (ROE_PPM(1) - ROE_IAM(1)) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, (ROE_IAM(1) - ROE_IAM(1)) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, (ROE_01(:,1) - ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_12(:,1) - ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_23(:,1) - ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_34(:,1) - ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_45(:,1) - ROE_IAM(1)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_56(:,1) - ROE_IAM(1)) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k', {'\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV', '\deltaV'})
hold off
grid on
legend('Start', 'Target', 'Maneuvering')
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(t_01(1) / T_chief, (ROE_PPM(2) - ROE_IAM(2)) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, (ROE_IAM(2) - ROE_IAM(2)) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, (ROE_01(:,2) - ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_12(:,2) - ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_23(:,2) - ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_34(:,2) - ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_45(:,2) - ROE_IAM(2)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_56(:,2) - ROE_IAM(2)) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
hold on
plot(t_01(1) / T_chief, (ROE_PPM(3) - ROE_IAM(3)) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, (ROE_IAM(3) - ROE_IAM(3)) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, (ROE_01(:,3) - ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_12(:,3) - ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_23(:,3) - ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_34(:,3) - ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_45(:,3) - ROE_IAM(3)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_56(:,3) - ROE_IAM(3)) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(t_01(1) / T_chief, (ROE_PPM(4) - ROE_IAM(4)) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, (ROE_IAM(4) - ROE_IAM(4)) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, (ROE_01(:,4) - ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_12(:,4) - ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_23(:,4) - ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_34(:,4) - ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_45(:,4) - ROE_IAM(4)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_56(:,4) - ROE_IAM(4)) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(t_01(1) / T_chief, (ROE_PPM(5) - ROE_IAM(5)) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, (ROE_IAM(5) - ROE_IAM(5)) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, (ROE_01(:,5) - ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_12(:,5) - ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_23(:,5) - ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_34(:,5) - ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_45(:,5) - ROE_IAM(5)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_56(:,5) - ROE_IAM(5)) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(t_01(1) / T_chief, (ROE_PPM(6) - ROE_IAM(6)) * a_chief * 1e3, 'ro')
plot(t_56(end) / T_chief, (ROE_IAM(6) - ROE_IAM(6)) * a_chief * 1e3, 'rx')
plot(t_01 / T_chief, (ROE_01(:,6) - ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_12 / T_chief, (ROE_12(:,6) - ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_23 / T_chief, (ROE_23(:,6) - ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_34 / T_chief, (ROE_34(:,6) - ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_45 / T_chief, (ROE_45(:,6) - ROE_IAM(6)) * a_chief * 1e3, 'b-')
plot(t_56 / T_chief, (ROE_56(:,6) - ROE_IAM(6)) * a_chief * 1e3, 'b-')
xline([t_01(1)/T_chief t_01(end)/T_chief t_12(end)/T_chief t_23(end)/T_chief t_34(end)/T_chief t_45(end)/T_chief], '--k')
hold off
grid on
ylabel('\deltaa [m]')
xlabel('Orbital Periods')


figure
subplot(3,1,1)
hold on
plot(ROE_PPM(2) * a_chief * 1e3, ROE_PPM(1) * a_chief * 1e3, 'ro', 'LineWidth', 2)
plot(ROE_IAM(2) * a_chief * 1e3, ROE_IAM(1) * a_chief * 1e3, 'rx', 'LineWidth', 2)
plot(ROE_01(:, 2) * a_chief * 1e3, ROE_01(:, 1) * a_chief * 1e3, 'b-')
plot(ROE_12(:, 2) * a_chief * 1e3, ROE_12(:, 1) * a_chief * 1e3, 'b-')
plot(ROE_23(:, 2) * a_chief * 1e3, ROE_23(:, 1) * a_chief * 1e3, 'b-')
plot(ROE_34(:, 2) * a_chief * 1e3, ROE_34(:, 1) * a_chief * 1e3, 'b')
plot(ROE_45(:, 2) * a_chief * 1e3, ROE_45(:, 1) * a_chief * 1e3, 'b')
plot(ROE_56(:, 2) * a_chief * 1e3, ROE_56(:, 1) * a_chief * 1e3, 'b')
hold off
grid on
ylabel('\deltaa [m]')
xlabel('\delta\lambda [m]')
legend('Start', 'Target', 'Maneuvering')

subplot(3,1,2)
hold on
plot(ROE_PPM(3) * a_chief * 1e3, ROE_PPM(4) * a_chief * 1e3, 'ro', 'LineWidth', 2)
plot(ROE_IAM(3) * a_chief * 1e3, ROE_IAM(4) * a_chief * 1e3, 'rx', 'LineWidth', 2)
scatter(ROE_01(:, 3) * a_chief * 1e3, ROE_01(:, 4) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_12(:, 3) * a_chief * 1e3, ROE_12(:, 4) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_23(:, 3) * a_chief * 1e3, ROE_23(:, 4) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_34(:, 3) * a_chief * 1e3, ROE_34(:, 4) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_45(:, 3) * a_chief * 1e3, ROE_45(:, 4) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_56(:, 3) * a_chief * 1e3, ROE_56(:, 4) * a_chief * 1e3, 'b', 'filled')
hold off
grid on
ylabel('\deltae_x [m]')
xlabel('\deltae_y [m]')

subplot(3,1,3)
hold on
plot(ROE_PPM(5) * a_chief * 1e3, ROE_PPM(6) * a_chief * 1e3, 'ro', 'LineWidth', 2)
plot(ROE_IAM(5) * a_chief * 1e3, ROE_IAM(6) * a_chief * 1e3, 'rx', 'LineWidth', 2)
scatter(ROE_01(:, 5) * a_chief * 1e3, ROE_01(:, 6) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_12(:, 5) * a_chief * 1e3, ROE_12(:, 6) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_23(:, 5) * a_chief * 1e3, ROE_23(:, 6) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_34(:, 5) * a_chief * 1e3, ROE_34(:, 6) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_45(:, 5) * a_chief * 1e3, ROE_45(:, 6) * a_chief * 1e3, 'b', 'filled')
scatter(ROE_56(:, 5) * a_chief * 1e3, ROE_56(:, 6) * a_chief * 1e3, 'b', 'filled')
hold off
grid on
ylabel('\deltai_x [m]')
xlabel('\deltai_y [m]')


%% Functions

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

function B = control_input(a, n, e, nu, omega, i)
    % for the usual quasi-non-singular ROE
    B = zeros(6,3);
    eta = sqrt(1-e^2);

    B(1,1) = 2 * e * sin(nu) / eta;
    B(1,2) = 2 * (1 + e * cos(nu)) / eta;

    B(2,1) = eta * (eta - 1) * cos(nu) / e - 2 * eta^2 / (1 + e * cos(nu));
    B(2,2) = eta * (1 - eta) * (2 + e * cos(nu)) * sin(nu) / (e * (1 + e * cos(nu)));

    B(3,1) = eta * sin(omega + nu);
    B(3,2) = eta * ((2 + e * cos(nu)) * cos(omega + nu) + e * cos(omega)) / (1 + e * cos(nu));
    B(3,3) = eta * e * sin(omega) * sin(omega + nu) / (tan(i) * (1 + e * cos(nu)));

    B(4,1) = -eta * cos(omega + nu);
    B(4,2) = eta * ((2 + e * cos(nu)) * sin(omega + nu) + e * sin(omega)) / (1 + e * cos(nu));
    B(4,3) = -eta * e * cos(omega) * sin(omega + nu) / (tan(i) * (1 + e * cos(nu)));
    
    B(5,3) = eta * cos(omega + nu) / (1 + e * cos(nu));

    B(6,3) = eta * sin(omega + nu) / (1 + e * cos(nu));

    B = B ./ a ./ n;
end

function statedot = GVE(t, state, mu)
    % for 2 satellites
    statedot = zeros(size(state));
    statedot(6) = sqrt(mu / state(1)^3);
    statedot(12) = sqrt(mu / state(7)^3);
    % a, e, i, omega, RAAN, M 
end

function delta_OE = transform(oe, delta_v, mu)
    delta_OE = zeros(1,6);
    % a, e, i, omega, RAAN, M
    nu = mean2true(oe(6), oe(2), 1e-12);
    n = sqrt(mu / oe(1)^3);
    r = oe(1) * (1 - oe(2)^2) / (1 + oe(2) * cos(nu));

    delta_OE(1) = 2 * oe(2) * sin(nu) / (n * sqrt(1 - oe(2)^2)) * delta_v(1) + 2 * oe(1) * sqrt(1 - oe(2)^2) / (n * r) * delta_v(2);
    delta_OE(2) = sqrt(1 - oe(2)^2) * sin(nu) / (n * oe(1)) * delta_v(1) + sqrt(1 - oe(2)^2) / (n * oe(1)^2 * oe(2)) * (oe(1)^2 * (1 - oe(2)^2) / r - r) * delta_v(2);
    delta_OE(3) = r * cos(oe(4) + nu) / (n * oe(1)^2 * sqrt(1 - oe(2)^2)) * delta_v(3);
    delta_OE(4) = -sqrt(1 - oe(2)^2) * cos(nu) / (n * oe(1) * oe(2)) * delta_v(1) + sqrt(1 - oe(2)^2) * (2 + oe(2) * cos(nu)) / (n * oe(1) * oe(2) * (1 + oe(2) * cos(nu))) * sin(nu) * delta_v(2) - r * cot(oe(3)) * sin(oe(4) + nu) / (n * oe(1)^2 * sqrt(1 - oe(2)^2)) * delta_v(3);
    delta_OE(5) = r * sin(oe(4) + nu) / (n * oe(1)^2 * sqrt(1 - oe(2)^2) * sin(oe(3))) * delta_v(3);
    delta_OE(6) = - 1 / (n * oe(1)) * (2 * r / oe(1) - (1 - oe(2)^2) * cos(nu) / oe(2)) * delta_v(1) - (1 - oe(2)^2) / (n * oe(1) * oe(2)) * (1 + r / (oe(1) * (1 - oe(2)^2))) * sin(nu) * delta_v(2);
end