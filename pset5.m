%% Problem 1

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
R_E = 6378; % km
J2 = 0.108263e-2;

% Defining the ROE for operational modes

% For the perigee pass
a_chief = 36943; % km
e_chief = 0.0001;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
M_chief = n_chief * (14 * 3600 + 49 * 60);
nu_chief = mean2true(M_chief, e_chief);
oe_chief = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

delta_a = 0;
delta_lambda = 0 / a_chief;
delta_ex = 800e-3 / a_chief;
delta_ey = 200e-3 / a_chief;
delta_ix = 750e-3 / a_chief;
delta_iy = 200e-3 / a_chief;

a_deputy = a_chief + delta_a * a_chief;
inc_deputy = inc_chief + delta_ix;
RAAN_deputy = RAAN_chief + delta_iy / sin(inc_chief);
e_deputy = sqrt((delta_ex + e_chief * cos(omega_chief))^2 + (delta_ey + e_chief * sin(omega_chief))^2);
omega_deputy = atan2(delta_ey + e_chief * sin(omega_chief), delta_ex + e_chief * cos(omega_chief));
M_deputy = delta_lambda + M_chief + omega_chief - (RAAN_deputy - RAAN_chief) * cos(inc_chief) - omega_deputy;
nu_deputy = mean2true(M_deputy, e_deputy);
oe_deputy = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

ROE = OE2ROE(oe_chief, oe_deputy) * a_chief * 1e3;

[pos_chief, vel_chief] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief, mu);
[pos_deputy, vel_deputy] = OE2ECI(a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, nu_deputy, mu);
rho_ECI = norm(pos_deputy - pos_chief);
norm(pos_chief);

% For the scientific experiment
M_chief = n_chief * (6 * 3600 + 49 * 60);
nu_chief = mean2true(M_chief, e_chief);
oe_chief = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

delta_a = 0;
delta_lambda = 0 / a_chief;
delta_ex = 59e-3 / a_chief;
delta_ey = 30e-3 / a_chief;
delta_ix = 80e-3 / a_chief;
delta_iy = 30e-3 / a_chief;

a_deputy = a_chief + delta_a * a_chief;
inc_deputy = inc_chief + delta_ix;
RAAN_deputy = RAAN_chief + delta_iy / sin(inc_chief);
e_deputy = sqrt((delta_ex + e_chief * cos(omega_chief))^2 + (delta_ey + e_chief * sin(omega_chief))^2);
omega_deputy = atan2(delta_ey + e_chief * sin(omega_chief), delta_ex + e_chief * cos(omega_chief));
M_deputy = delta_lambda + M_chief + omega_chief - (RAAN_deputy - RAAN_chief) * cos(inc_chief) - omega_deputy;
nu_deputy = mean2true(M_deputy, e_deputy);
oe_deputy = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

ROE = OE2ROE(oe_chief, oe_deputy) * a_chief * 1e3;

[pos_chief, vel_chief] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief, mu);
[pos_deputy, vel_deputy] = OE2ECI(a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, nu_deputy, mu);
rho_ECI = norm(pos_deputy - pos_chief);
norm(pos_chief);

% For the formation flying experiment
M_chief = n_chief * (6 * 3600 + 49 * 60);
nu_chief = mean2true(M_chief, e_chief);
oe_chief = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];

delta_a = 0;
delta_lambda = 0 / a_chief;
delta_ex = 50e-3 / a_chief;
delta_ey = 89e-3 / a_chief;
delta_ix = 100e-3 / a_chief;
delta_iy = 120e-3 / a_chief;

a_deputy = a_chief + delta_a * a_chief;
inc_deputy = inc_chief + delta_ix;
RAAN_deputy = RAAN_chief + delta_iy / sin(inc_chief);
e_deputy = sqrt((delta_ex + e_chief * cos(omega_chief))^2 + (delta_ey + e_chief * sin(omega_chief))^2);
omega_deputy = atan2(delta_ey + e_chief * sin(omega_chief), delta_ex + e_chief * cos(omega_chief));
M_deputy = delta_lambda + M_chief + omega_chief - (RAAN_deputy - RAAN_chief) * cos(inc_chief) - omega_deputy;
nu_deputy = mean2true(M_deputy, e_deputy);
oe_deputy = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];

ROE = OE2ROE(oe_chief, oe_deputy) * a_chief * 1e3

[pos_chief, vel_chief] = OE2ECI(a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, nu_chief, mu);
[pos_deputy, vel_deputy] = OE2ECI(a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, nu_deputy, mu);
rho_ECI = norm(pos_deputy - pos_chief)

% Funtions

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end