%% Problem 2

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
R_E = 6378; % km
J2 = 0.108263e-2;

options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100); % options for numerical integration

% Comparing integrating the position/velocity in ECI vs ROE over 100
% orbits

% Initial conditions = the ones at the start of the perigee pass mode
a_chief = 36943; % km
e_chief = 0.8111;
% e_chief = 1e-5;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
T_chief = 2 * pi / n_chief;
M_chief = n_chief * (14 * 3600 + 49 * 60); % perigee pass mode starts at t=14h49
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

ROE_PPM = OE2ROE(oe_chief_PPM, oe_deputy_PPM);

[pos_chief, vel_chief] = OE2ECI(oe_chief_PPM(1), oe_chief_PPM(2), oe_chief_PPM(3), oe_chief_PPM(4), oe_chief_PPM(5), nu_chief_PPM, mu);
[pos_deputy, vel_deputy] = OE2ECI(oe_deputy_PPM(1), oe_deputy_PPM(2), oe_deputy_PPM(3), oe_deputy_PPM(4), oe_deputy_PPM(5), nu_deputy_PPM, mu);
ic_ECI = [pos_chief', vel_chief', pos_deputy', vel_deputy'];

% Numerical integration of the different sets of equations over 100 orbits

N = 100000;
tspan = linspace(0, 100 * T_chief, N);

[t_ECI, y_ECI] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, ic_ECI, options);
% [t_ECI, y_ECI] = ode89(@(t, state) FODE_2sats_J2(t, state, mu, J2, R_E), tspan, ic_ECI, options);
[t_ROE, y_ROE] = ode89(@(t, state) ROE_noJ2(t, state, mu, a_chief, a_deputy), tspan, ROE_PPM, options);

% Computing the ROE at each time step from the propagation of the
% position/velocity state in ECI

ROE_history = zeros(N, 6);
ROE_history_STM = zeros(N, 6);

for j=1:N
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j,1:6), mu);
    oe_chief_temp = [a, e, i, omega, RAAN, M];
%     nu_chief_temp = mean2true(M, e);
    [a, e, i, omega, RAAN, M] = ECI2OE_M(y_ECI(j,7:12), mu);
    oe_deputy_temp = [a, e, i, omega, RAAN, M];
    ROE_cur = OE2ROE(oe_chief_temp, oe_deputy_temp);
    ROE_history(j, :) = ROE_cur;
    A = ROE_STM(tspan(j), sqrt(mu / oe_chief_temp(1)^3));
    ROE_cur_STM = A * ROE_PPM';
    ROE_history_STM(j,:) = ROE_cur_STM;
end

figure
subplot(3,2,1)
plot(tspan / T_chief, (ROE_history(:,1) - y_ROE(:,1)) * a_chief * 1e3)
ylabel('\Delta\deltaa [m]')

subplot(3,2,3)
plot(tspan / T_chief, (ROE_history(:,2) - y_ROE(:,2)) * a_chief * 1e3)
ylabel('\Delta\delta\lambda [m]')

subplot(3,2,5)
plot(tspan / T_chief, (ROE_history(:,3) - y_ROE(:,3)) * a_chief * 1e3)
ylabel('\Delta\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
plot(tspan / T_chief, (ROE_history(:,4) - y_ROE(:,4)) * a_chief * 1e3)
ylabel('\Delta\deltae_y [m]')

subplot(3,2,4)
plot(tspan / T_chief, (ROE_history(:,5) - y_ROE(:,5)) * a_chief * 1e3)
ylabel('\Delta\deltai_x [m]')

subplot(3,2,6)
plot(tspan / T_chief, (ROE_history(:,6) - y_ROE(:,6)) * a_chief * 1e3)
ylabel('\Delta\deltai_y [m]')
xlabel('Orbital Periods')


figure
subplot(3,2,1)
plot(tspan / T_chief, (ROE_history(:,1) - ROE_history_STM(:,1)) * a_chief * 1e3)
ylabel('\Delta\deltaa [m]')

subplot(3,2,3)
plot(tspan / T_chief, (ROE_history(:,2) - ROE_history_STM(:,2)) * a_chief * 1e3)
ylabel('\Delta\delta\lambda [m]')

subplot(3,2,5)
plot(tspan / T_chief, (ROE_history(:,3) - ROE_history_STM(:,3)) * a_chief * 1e3)
ylabel('\Delta\deltae_x [m]')
xlabel('Orbital Periods')

subplot(3,2,2)
plot(tspan / T_chief, (ROE_history(:,4) - ROE_history_STM(:,4)) * a_chief * 1e3)
ylabel('\Delta\deltae_y [m]')

subplot(3,2,4)
plot(tspan / T_chief, (ROE_history(:,5) - ROE_history_STM(:,5)) * a_chief * 1e3)
ylabel('\Delta\deltai_x [m]')

subplot(3,2,6)
plot(tspan / T_chief, (ROE_history(:,6) - ROE_history_STM(:,6)) * a_chief * 1e3)
ylabel('\Delta\deltai_y [m]')
xlabel('Orbital Periods')


%% Functions

function statedot = ROE_noJ2(t, state, mu, a_c, a_d)
    statedot = zeros(size(state));
    statedot(2) = sqrt(mu) * (1/sqrt(a_d^3) - 1/sqrt(a_c^3));
end

function statedot = ROE_J2(t, state, mu, a_c, a_d, i_c, i_d, e_c, e_d, omega_c, omega_d)
    statedot = zeros(size(state));
    
    eta_c = sqrt(1-e_c^2);
    eta_d = sqrt(1 - e_d^2);
    k_c = 3 * J2 * R_E^2 * sqrt(mu) / (4 * a_c^(7/2) * eta_c^4);
    k_d = 3 * J2 * R_E^2 * sqrt(mu) / (4 * a_d^(7/2) * eta_d^4);

    der_d = k_d * [0, (eta_d *(3*cos(i_d)^2 - 1) + (5*cos(i_d)^2 - 1) - 2*cos(i_d)*cos(i_c)), ...
        -e_d * (5*cos(i_d)^2 - 1)*sin(omega_d), e_d*(5*cos(i_d)^2 - 1)*cos(omega_d), 0, -2*cos(i_d)*sin(i_c)];
    der_c = k_c * [0, (eta_c *(3*cos(i_c)^2 - 1) + (5*cos(i_c)^2 - 1) - 2*cos(i_c)^2), ...
        -e_c * (5*cos(i_c)^2 - 1)*sin(omega_c), e_c*(5*cos(i_c)^2 - 1)*cos(omega_chief), 0, -2*cos(i_c)*sin(i_c)];
    statedot(:) = der_d - der_c + sqrt(mu) * [0, 1/sqrt(a_d^3) - 1/sqrt(a_c^3), 0, 0, 0, 0];

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

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end

function A = ROE_STM(t, n)
    A = eye(6);
    A(2,1) = -1.5 * n * t;
end