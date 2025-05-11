%% Problem 1

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
% R_E = 6378; % km
% J2 = 0.108263e-2;
T_max = 10e-6 / 340; % maximum along track acceleration
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

% Setting up the simulation for reconfiguration using continuous control

[pos_chief, vel_chief] = OE2ECI(oe_chief_PPM(1), oe_chief_PPM(2), oe_chief_PPM(3), oe_chief_PPM(4), oe_chief_PPM(5), nu_chief_PPM, mu);
[pos_deputy, vel_deputy] = OE2ECI(oe_deputy_PPM(1), oe_deputy_PPM(2), oe_deputy_PPM(3), oe_deputy_PPM(4), oe_deputy_PPM(5), nu_deputy_PPM, mu);

N = 40000;
reconfiguration_tspan = linspace(0, 40 * T_chief, N); % to first try and debug code before doing it for longer to see convergence
state_chief_history = zeros(N,6);
state_chief_history(1,:) = [pos_chief', vel_chief'];

state_deputy_history = zeros(N,6);
state_deputy_history(1,:) = [pos_deputy', vel_deputy'];

modified_ROE_history = zeros(N, 6);
u_history = zeros(N, 2);
control_tracking_error_history = zeros(N,6);

% Do a first try where you don't control delta_lambda but everything else,
% then add the control of delta_lambda, then add the reference governor if
% time (but very difficult)

for j=1:N-1
%     j
    % first compute the control law
    % Comptue the current modified ROE
    [a, e, i, omega, RAAN, M] = ECI2OE_M(state_chief_history(j,:), mu);
    oe_chief_temp = [a, e, i, omega, RAAN, M];
    nu_chief_temp = mean2true(M, e);
    [a, e, i, omega, RAAN, M] = ECI2OE_M(state_deputy_history(j,:), mu);
    oe_deputy_temp = [a, e, i, omega, RAAN, M];
    nu_deputy_temp = mean2true(M, e);
    modified_ROE_cur = eccentric_ROE(oe_chief_temp, oe_deputy_temp);
    ROE_cur = OE2ROE(oe_chief_temp, oe_deputy_temp);
    modified_ROE_history(j,:) = modified_ROE_cur;
    delta_alpha = modified_ROE_cur - modified_ROE_IAM;
    control_tracking_error_history(j,:) = delta_alpha;

    % For now, not controlling delta_lambda so applied reference =
    % reference
    modified_ROE_applied = modified_ROE_IAM;

    % Trying to control delta_lambda
    if abs(delta_alpha(2)) > 1e-6
%         j
        % From Steindorf
%         delta_v_orbit = 3 / 4 * pi * T_max; % specifically when N = 4 in gain matrix
%         delta_a_tan = 2 / (oe_chief_temp(1) * sqrt(mu / oe_chief_temp(1)^3) * sqrt(1-oe_chief_temp(2)^2)) * (1+oe_chief_temp(2)*cos(nu_chief_temp)) * delta_v_orbit / 2;
%         delta_a_des = abs(delta_a_tan) / 2;
%         delta_lambda_ref = sign(delta_alpha(2)) * min([abs(delta_alpha(2)) / (100 * T_chief), 3 / 2 * sqrt(mu / oe_chief_temp(1)^3) * delta_a_des]);
%         delta_a_applied = -2 / 3 * delta_lambda_ref / sqrt(mu / oe_chief_temp(1^3));
%         modified_ROE_applied(1) = delta_a_applied;

        % From Lippe
        delta_a_applied = sign(delta_alpha(2)) * min([abs(delta_alpha(2)) / (500e-3 / oe_chief_temp(1)), 200e-3 / oe_chief_temp(1)]);
        modified_ROE_applied(1) = delta_a_applied;
    end

    delta_alpha_applied = modified_ROE_cur - modified_ROE_applied;
%     delta_alpha_applied = delta_alpha;

    % Don't forget to get rid off the delta lambda before doing the rest
    B = modified_control_matrix(oe_chief_temp(1), sqrt(mu / oe_chief_temp(1)^3), oe_chief_temp(2), nu_chief_temp, oe_chief_temp(4));
    A = zeros(5);

    % Computing the optimal location of the out-of-plane maneuver
    nu_oop_temp = atan2(delta_alpha(6), delta_alpha(5));
    if cos(nu_oop_temp) <= cos(nu_oop_temp + pi)
        nu_oop = nu_oop_temp;
    else
        nu_oop = nu_oop_temp + pi;
    end
    delta_nu_oop = pi; % it is what is done for near circular but I can adapt it

    % Computing the optimal location for the in-plane maneuver
    delta_ex_tild = delta_alpha_applied(3);
    delta_ey_tild = oe_chief_temp(2) * delta_alpha_applied(4);
%     delta_ey_tild = delta_alpha_applied(4);
    if delta_ex_tild == 0
        nu_ip = acos((sqrt(1-oe_chief_temp(2)^2) - 1) / oe_chief_temp(2));
    elseif delta_ey_tild == 0
        nu_ip = 0;
    else
        % make the function that will return the in-plane optimal maneuver
        % time
        [nu_ip1, nu_ip2] = inplane_maneuver(delta_ex_tild, delta_ey_tild, oe_chief_temp(2), sqrt(mu / oe_chief_temp(1)^3));
%         if abs(nu_ip1 - nu_deputy_temp) <= abs(nu_ip2 - nu_deputy_temp)
%             nu_ip = nu_ip1;
%         else
%             nu_ip = nu_ip2;
%         end
%         nu_ip = atan2(delta_ey_tild, delta_ex_tild);
        nu_ip = nu_ip1;
    end
%     nu_ip, nu_chief_temp
    delta_nu_ip = oe_chief_temp(1) * sqrt(mu / oe_chief_temp(1)^3) * sqrt(1-oe_chief_temp(2)^2) / (2 * T_max) * sqrt(mu / (oe_chief_temp(1)^3 * (1-oe_chief_temp(2)^2)^3)) * 1 / (1 + oe_chief_temp(2) * cos(nu_ip))^3 * abs(delta_alpha_applied(1));

%     P = control_gain(3000, nu_deputy_temp, nu_ip, nu_oop, delta_nu_ip, delta_nu_oop);
    P = control_gain2(3000, 4, nu_deputy_temp, nu_ip, nu_oop);

    % Finally computing the control
%     u = - (B' * B)^(-1) * B' * (A * [modifi ed_ROE_cur(1); modified_ROE_cur(3:end)'] + P * [delta_alpha_applied(1); delta_alpha_applied(3:end)']);
    u = - pinv(B) * (A * [modified_ROE_cur(1); modified_ROE_cur(3:end)'] + P * [delta_alpha_applied(1); delta_alpha_applied(3:end)']);
    % should be 2D because no radial maneuvers 
    u_history(j,:) = u;

    % Rotate the delta v from RTN to ECI (using deputy's RTN but look at Steindorf's work)
    % RTN unit vectors
    R = state_deputy_history(j, 1:3) / norm(state_deputy_history(j, 1:3));
    h = cross(state_deputy_history(j, 1:3), state_deputy_history(j, 4:6));
    N = h / norm(h);
    T = cross(N, R);
    
    % Rotation matrix from ECI to RTN
    rotation = [R; T; N];
    delta_v = rotation' * [0; u];

    % Apply the delta v to the deputy
%     new_state_deputy = state_deputy_history(j,:) + [0 0 0 delta_v'];

    % Numerically integrate until the next time step and store it in a
    % matrix to be able to access it in the next loop for control law
    tspan = [reconfiguration_tspan(j), reconfiguration_tspan(j+1)];
    ic = [state_chief_history(j, :), state_deputy_history(j,:)];
    [t, y] = ode89(@(t, state) FODE_2sats(t, state, mu, delta_v), tspan, ic, options);
    state_chief_history(j+1,:) = y(end, 1:6);
    state_deputy_history(j+1, :) = y(end, 7:12);
end

[a, e, i, omega, RAAN, M] = ECI2OE_M(state_chief_history(end,:), mu);
oe_chief_temp = [a, e, i, omega, RAAN, M];
% nu_chief_temp = mean2true(M, e);
[a, e, i, omega, RAAN, M] = ECI2OE_M(state_deputy_history(end,:), mu);
oe_deputy_temp = [a, e, i, omega, RAAN, M];
% nu_deputy_temp = mean2true(M, e);
modified_ROE_cur = eccentric_ROE(oe_chief_temp, oe_deputy_temp);
modified_ROE_history(end,:) = modified_ROE_cur;
control_tracking_error_history(end,:) = modified_ROE_cur - modified_ROE_IAM;


figure
subplot(3,2,1)
hold on
plot(0, modified_ROE_PPM(1) * a_chief * 1e3, 'ro')
plot(reconfiguration_tspan / T_chief, modified_ROE_history(:,1) * a_chief * 1e3, 'b-')
plot(reconfiguration_tspan(end) / T_chief, modified_ROE_IAM(1) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltaa [m]')

subplot(3,2,3)
hold on
plot(0, modified_ROE_PPM(2) * a_chief * 1e3, 'ro')
plot(reconfiguration_tspan / T_chief, modified_ROE_history(:,2) * a_chief * 1e3, 'b-')
plot(reconfiguration_tspan(end) / T_chief, modified_ROE_IAM(2) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\delta\lambda_e [m]')

subplot(3,2,5)
hold on
plot(0, modified_ROE_PPM(3) * a_chief * 1e3, 'ro')
plot(reconfiguration_tspan / T_chief, modified_ROE_history(:,3) * a_chief * 1e3, 'b-')
plot(reconfiguration_tspan(end) / T_chief, modified_ROE_IAM(3) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltae_x [m]')

subplot(3,2,2)
hold on
plot(0, modified_ROE_PPM(4) * a_chief * 1e3, 'ro')
plot(reconfiguration_tspan / T_chief, modified_ROE_history(:,4) * a_chief * 1e3, 'b-')
plot(reconfiguration_tspan(end) / T_chief, modified_ROE_IAM(4) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
hold on
plot(0, modified_ROE_PPM(5) * a_chief * 1e3, 'ro')
plot(reconfiguration_tspan / T_chief, modified_ROE_history(:,5) * a_chief * 1e3, 'b-')
plot(reconfiguration_tspan(end) / T_chief, modified_ROE_IAM(5) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
hold on
plot(0, modified_ROE_PPM(6) * a_chief * 1e3, 'ro')
plot(reconfiguration_tspan / T_chief, modified_ROE_history(:,6) * a_chief * 1e3, 'b-')
plot(reconfiguration_tspan(end) / T_chief, modified_ROE_IAM(6) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel('\deltai_y [m]')

figure
subplot(3,2,1)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,1) * a_chief * 1e3)
ylabel('\deltaa [m]')

subplot(3,2,3)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,2) * a_chief * 1e3)
ylabel('\delta\lambda [m]')

subplot(3,2,5)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,3) * a_chief * 1e3)
ylabel('\deltae_x [m]')

subplot(3,2,2)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,4) * a_chief * 1e3)
ylabel('\deltae_y [m]')

subplot(3,2,4)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,5) * a_chief * 1e3)
ylabel('\deltai_x [m]')

subplot(3,2,6)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,6) * a_chief * 1e3)
ylabel('\deltai_y [m]')

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

function statedot = FODE_2sats(t, state, mu, control)
    statedot = zeros(size(state));
    statedot(1:3) = state(4:6);
    statedot(7:9) = state(10:12);

    r0 = state(1:3);
    r1 = state(7:9);

    statedot(4:6) = - mu * r0 / norm(r0)^3;
    statedot(10:12) = - mu * r1 / norm(r1)^3 + control;
end

function B = modified_control_matrix(a, n, e, nu, omega)
%     B = zeros(6, 3);
    B = zeros(5,2);
    eta = sqrt(1-e^2);

%     B(1,1) = 2 * e * sin(nu) / eta;
%     B(1,2) = 2 * (1+e*cos(nu)) / eta;
    B(1,1) = 2 * (1+e*cos(nu)) / eta;

%     B(2,1) = - 2 * eta^2 / (1+e*cos(nu));

%     B(3,1) = eta * sin(nu);
%     B(3,2) = eta * (e + cos(nu)*(2+e*cos(nu))) / (1+e*cos(nu));
    B(2,1) = eta * (e + cos(nu)*(2+e*cos(nu))) / (1+e*cos(nu));
    
%     B(4,1) = - eta * cos(nu) / e;
%     B(4,2) = (eta / e) * sin(nu) * (2+e*cos(nu)) / (1+e*cos(nu));
    B(3,1) = (eta / e) * sin(nu) * (2+e*cos(nu)) / (1+e*cos(nu));

%     B(5,3) = eta * cos(omega + nu) / (1+e*cos(nu));
    B(4,2) = eta * cos(omega + nu) / (1+e*cos(nu));

%     B(6,3) = eta * sin(omega + nu) / (1+e*cos(nu));
    B(5,2) = eta * sin(omega + nu) / (1+e*cos(nu));

    B = B ./ (a .* n);
end

function A = STM(n, delta_t)
    A = eye(6);
    A(2,1) = - 1.5 * n * delta_t;
end

function P = control_gain(k, nu, nu_ip, nu_oop, delta_nu_ip, delta_nu_oop)
%     P = zeros(6,6);
    P = zeros(5,5);
    if abs(nu - nu_ip) < delta_nu_ip/2
        P(1,1) = 1;
        P(2,2) = 1;
        P(3,3) = 1;
%         P(4,4) = 1;
    end
    if abs(nu - nu_oop) < delta_nu_oop/2
        P(4,4) = 1;
        P(5,5) = 1;
%         P(6,6) = 1;
    end
    P = (1 ./ k) .* P;
end

function P = control_gain2(k, N, nu, nu_ip, nu_oop)
    P = zeros(5);
    P(1,1) = cos(nu - nu_ip)^N;
    P(2,2) = cos(nu - nu_ip)^N;
    P(3,3) = cos(nu - nu_ip)^N;
    P(4,4) = cos(nu - nu_oop)^N;
    P(5,5) = cos(nu - nu_oop)^N;
    P = P ./ k;
end

function [nu_ip1, nu_ip2] = inplane_maneuver(delta_ex_tild, delta_ey_tild, e, n)
    nu_dis = wrapToPi(pi + acos(e));
    nu_re = wrapToPi(pi - acos(e));
    phi = atan2(delta_ey_tild, delta_ex_tild); % phase of the eccentricity vector of pseudostate

    if delta_ex_tild * delta_ey_tild > 0
        if ((nu_dis < phi) && (phi < wrapToPi(pi - nu_dis))) || ((wrapToPi(pi - nu_re) < phi) && (phi < nu_re))
            range1 = [nu_dis, 2*pi];
            range2 = [0, nu_re];
        else
            range1 = [pi, nu_dis];
            range2 = [0, nu_re];
        end
    elseif delta_ex_tild * delta_ey_tild < 0
        if ((nu_dis < phi) && (phi < wrapToPi(pi - nu_dis))) || ((wrapToPi(pi - nu_re) < phi) && (phi < nu_re))
            range1 = [0, nu_re];
            range2 = [nu_dis, 2*pi];
        else
            range1 = [nu_re, pi];
            range2 = [nu_dis, 2*pi];
        end
    end
    
    nu_ip1 = wrapToPi(NR((range1(1) + range1(2)) / 2, e, n, delta_ex_tild, delta_ey_tild));
    nu_ip2 = wrapToPi(NR((range2(1) + range2(2)) / 2, e, n, delta_ex_tild, delta_ey_tild));
end

function x = NR(init_guess, e, n, delta_ex_des, delta_ey_des)
%     x_prev = inf;
    x = init_guess;
    i = 0;
    f_val = f(x, e, n, delta_ex_des, delta_ey_des);

    while abs(f_val) > 1e-12 && i < 100
%         f(x, e, n, delta_ex_des, delta_ey_des)
%         f_der(x, e, n, delta_ex_des, delta_ey_des)
        new_x = x - f(x, e, n, delta_ex_des, delta_ey_des) / f_der(x, e, n, delta_ex_des, delta_ey_des);
%         x_prev = x;
        x = new_x;
        f_val = f(x, e, n, delta_ex_des, delta_ey_des);
        i = i+1;
    end
%     i
end

function x = f(nu, e, n, delta_ex_des, delta_ey_des)
    delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des);
    delta_ex_tild = sqrt(1-e^2) / n * (sin(nu) * sqrt(1 - delta_v^2) + ((2+e*cos(nu)) * cos(nu) + e) / (1+e*cos(nu)) * delta_v);
    delta_ey_tild = sqrt(1-e^2) / n * (-cos(nu) * sqrt(1 - delta_v^2) + ((2+e*cos(nu)) * sin(nu)) / (1+e*cos(nu)) * delta_v);
    x = delta_ey_tild - delta_ey_des / delta_ex_des * delta_ex_tild;
end

function delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des)
    if sign(delta_ex_des) == sign(delta_ey_des)
        delta_v = - sqrt(0.5 - f1(nu, e) / (2 * sqrt(4 + f1(nu, e)^2)));
    else
        delta_v = sqrt(0.5 + f1(nu, e) / (2 * sqrt(4 + f1(nu, e)^2)));
    end
end

function num = f1(nu, e)
    num = f2(nu, e) / ((1  + e * cos(nu)) * e * sin(nu));
end

function x = f2(nu, e)
    x = 2 * e^2 * cos(nu)^2 + 6 * e * cos(nu) + e^2 + 3;
end

function x = f_der(nu, e, n, delta_ex_des, delta_ey_des)
    der_delta_ex_tild_num = der_delta_ex_tild(nu, e, n, delta_ex_des, delta_ey_des);
    der_delta_ey_tild_num = der_delta_ey_tild(nu, e, n, delta_ex_des, delta_ey_des);
    x = der_delta_ey_tild_num - delta_ey_des / delta_ex_des * der_delta_ex_tild_num;
end

function x = der_delta_ex_tild(nu, e, n, delta_ex_des, delta_ey_des)
    delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des);
    der_delta_v = der_delta_vt(nu, e, delta_ex_des, delta_ey_des);
    x = sqrt(1-e^2) / n * (cos(nu) * sqrt(1-delta_v^2) - sin(nu) * der_delta_v * delta_v * (1-delta_v^2) ...
        + (delta_v * (-sin(nu) * (2+e*cos(nu)) - cos(nu) * sin(nu) * e) + der_delta_v * ((2+e*cos(nu)) * cos(nu)+e))*(1+e*cos(nu)) / (1+e*cos(nu))^2 ...
        + delta_v * ((2+e*cos(nu))*cos(nu)+e) * e * sin(nu) / (1+e*cos(nu))^2);
end

function x = der_delta_ey_tild(nu, e, n, delta_ex_des, delta_ey_des)
    delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des);
    der_delta_v = der_delta_vt(nu, e, delta_ex_des, delta_ey_des);
    x = sqrt(1-e^2) / n * (sin(nu)*sqrt(1-delta_v^2) + cos(nu)*der_delta_v*delta_v*(1-delta_v^2) ...
        + (der_delta_v*sin(nu)*(2+e*cos(nu)) +  delta_v*(cos(nu) * (2+e*cos(nu)) - sin(nu)^2*e)) * (1+e*cos(nu)) / (1+e*cos(nu))^2 ...
        + (2+e*cos(nu))*sin(nu)*delta_v*e*sin(nu) / (1+e*cos(nu))^2);
end

function x = der_delta_vt(nu, e, delta_ex_des, delta_ey_des)
    % need to add the case when the deltas are of opposite signs
    delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des);
    f1_num = f1(nu, e);
    der_f1_num = der_f1(nu, e);
    x = 0.5 * (-1 / delta_v) * (der_f1_num * 2 * (4+f1_num^2)^(0.5) / (4 * (4 +f1_num^2)) ...
        - f1_num^2 * der_f1_num * (4+f1_num^2)^(-0.5) / (4*(4+f1_num^2)));
end

function x = der_f1(nu, e)
    f2_num = f2(nu, e);
    der_f2_num = der_f2(nu, e);
    x = (der_f2_num * (1+e*cos(nu)) * e * sin(nu) - f2_num * (e*cos(nu)*(1+e*cos(nu)) - e^2*sin(nu)^2)) / ((1+e*cos(nu))^2 *e^2 * sin(nu)^2);
end

function x = der_f2(nu, e)
    x = -4*e^2*cos(nu)*sin(nu) - 6*e*sin(nu);
end

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapTo2Pi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = (oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3));
end

function oe_deputy = ROE2OE(oe_chief, ROE)
    % if time to write it later
end