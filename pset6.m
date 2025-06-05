%% Problem 1

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
% T_max = 10e-6 / 340; % maximum along track acceleration
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

N = 20000;
reconfiguration_tspan = linspace(0, 20 * T_chief, N);
dt = reconfiguration_tspan(2) - reconfiguration_tspan(1);
state_chief_history = zeros(N,6);
state_chief_history(1,:) = [pos_chief', vel_chief'];

state_deputy_history = zeros(N,6);
state_deputy_history(1,:) = [pos_deputy', vel_deputy'];

modified_ROE_history = zeros(N, 6);
ROE_history = zeros(N, 6);
u_history = zeros(N, 2);
total_dv = zeros(N,2);
control_tracking_error_history = zeros(N,6);

% Do a first try where you don't control delta_lambda but everything else,
% then add the control of delta_lambda, then add the reference governor if
% time (but very difficult)

for j=1:N-1
    j
    % first compute the control law
    % Compute the current modified ROE
    [a, e, i, omega, RAAN, M] = ECI2OE_M(state_chief_history(j,:), mu);
    oe_chief_temp = [a, e, i, omega, RAAN, M];
    nu_chief_temp = mean2true(M, e);
    [a, e, i, omega, RAAN, M] = ECI2OE_M(state_deputy_history(j,:), mu);
    oe_deputy_temp = [a, e, i, omega, RAAN, M];
    nu_deputy_temp = mean2true(M, e);
    modified_ROE_cur = eccentric_ROE(oe_chief_temp, oe_deputy_temp);
    ROE_cur = OE2ROE(oe_chief_temp, oe_deputy_temp);
    ROE_history(j,:) = ROE_cur;
    modified_ROE_history(j,:) = modified_ROE_cur;
    delta_alpha = modified_ROE_cur - modified_ROE_IAM;
%     delta_alpha = ROE_cur - ROE_IAM;
    control_tracking_error_history(j,:) = delta_alpha;

    % If not controlling delta_lambda, applied = goal
    modified_ROE_applied = modified_ROE_IAM;
%     modified_ROE_applied = ROE_IAM;

    % Controlling delta_lambda
    if abs(delta_alpha(2)) > 1e-6
        % From Lippe
        delta_a_applied = sign(delta_alpha(2)) * min([abs(delta_alpha(2)) / (500e-3 / oe_chief_temp(1)), 100e-3 / oe_chief_temp(1)]);
        modified_ROE_applied(1) = delta_a_applied;
    end

    delta_alpha_applied = modified_ROE_cur - modified_ROE_applied;

    B = modified_control_matrix(oe_chief_temp(1), sqrt(mu / oe_chief_temp(1)^3), oe_chief_temp(2), nu_chief_temp, oe_chief_temp(4));
%     B = control_matrix(oe_chief_temp(1), sqrt(mu / oe_chief_temp(1)^3), oe_chief_temp(2), nu_chief_temp, oe_chief_temp(4), oe_chief_temp(3));
    A = zeros(5);

    % Computing the optimal location of the out-of-plane maneuver
    nu_oop = wrapToPi(atan2(delta_alpha(6), delta_alpha(5)) - oe_chief_temp(4));

    % Computing the optimal location for the in-plane maneuver
    delta_ex_tild = delta_alpha_applied(3);
    delta_ey_tild = oe_chief_temp(2) * delta_alpha_applied(4);
    if delta_ex_tild == 0
        nu_ip = acos((sqrt(1-oe_chief_temp(2)^2) - 1) / oe_chief_temp(2));
    elseif delta_ey_tild == 0
        nu_ip = 0;
    else
        del_ex = oe_deputy_temp(2) * cos(oe_deputy_temp(4)) - oe_chief_temp(2) * cos(oe_chief_temp(4)) - (oe_deputy_IAM(2) * cos(oe_deputy_IAM(4)) - oe_chief_IAM(2) * cos(oe_chief_IAM(4)));
        del_ey = oe_deputy_temp(2) * sin(oe_deputy_temp(4)) - oe_chief_temp(2) * sin(oe_chief_temp(4)) - (oe_deputy_IAM(2) * sin(oe_deputy_IAM(4)) - oe_chief_IAM(2) * sin(oe_chief_IAM(4)));
        nu_ip = cast(wrapToPi(inplane_planning(oe_chief_temp(2), oe_chief_temp(4), oe_chief_temp(3), del_ex, del_ey, delta_alpha_applied(6)) - oe_chief_temp(4)), 'double');
    end

    P = control_gain(3000, 4, nu_deputy_temp, nu_ip, nu_oop); % try different values for N

    % Finally computing the control
    u = - pinv(B) * (A * [modified_ROE_cur(1); modified_ROE_cur(3:end)'] + P * [delta_alpha_applied(1); delta_alpha_applied(3:end)']);
    % should be 2D because no radial maneuvers 
    u_history(j+1,:) = u; % in km/s^2
    total_dv(j+1,1) = sum(abs(u_history(1:j, 1)) * dt * 1e3); % in m/s
    total_dv(j+1,2) = sum(abs(u_history(1:j, 2)) * dt * 1e3);

    % Rotate the delta v from RTN to ECI (using deputy's RTN)
    % RTN unit vectors
    R = state_deputy_history(j, 1:3) / norm(state_deputy_history(j, 1:3));
    h = cross(state_deputy_history(j, 1:3), state_deputy_history(j, 4:6));
    N = h / norm(h);
    T = cross(N, R);
    
    % Rotation matrix from ECI to RTN
    rotation = [R; T; N];
    delta_v = rotation' * [0; u];

    % Numerically integrate until the next time step and store it in a
    % matrix to be able to access it in the next loop for control law
    tspan = [reconfiguration_tspan(j), reconfiguration_tspan(j+1)];
    ic = [state_chief_history(j, :), state_deputy_history(j,:)];
    [t, y] = ode89(@(t, state) FODE_2sats(t, state, mu, delta_v), tspan, ic, options);
    state_chief_history(j+1,:) = y(end, 1:6);
    state_deputy_history(j+1, :) = y(end, 7:12);
end

% Computing the ROE after the last integration
[a, e, i, omega, RAAN, M] = ECI2OE_M(state_chief_history(end,:), mu);
oe_chief_temp = [a, e, i, omega, RAAN, M];
[a, e, i, omega, RAAN, M] = ECI2OE_M(state_deputy_history(end,:), mu);
oe_deputy_temp = [a, e, i, omega, RAAN, M];
modified_ROE_cur = eccentric_ROE(oe_chief_temp, oe_deputy_temp);
ROE_cur = OE2ROE(oe_chief_temp, oe_deputy_temp);
ROE_history(end,:) = ROE_cur;
modified_ROE_history(end,:) = modified_ROE_cur;
control_tracking_error_history(end,:) = modified_ROE_cur - modified_ROE_IAM;
% control_tracking_error_history(end,:) = ROE_cur - ROE_IAM;

%% Plotting

figure
subplot(3,2,1)
hold on
plot(0, modified_ROE_PPM(1) * a_chief * 1e3, 'ro')
plot(reconfiguration_tspan / T_chief, modified_ROE_history(:,1) * a_chief * 1e3, 'b-')
plot(reconfiguration_tspan(end) / T_chief, modified_ROE_IAM(1) * a_chief * 1e3, 'rx')
hold off
legend('Start', 'Maneuvering', 'Target')
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
ylabel("\deltae_x' [m]")
xlabel('Orbital Periods')

subplot(3,2,2)
hold on
plot(0, modified_ROE_PPM(4) * a_chief * 1e3, 'ro')
plot(reconfiguration_tspan / T_chief, modified_ROE_history(:,4) * a_chief * 1e3, 'b-')
plot(reconfiguration_tspan(end) / T_chief, modified_ROE_IAM(4) * a_chief * 1e3, 'rx')
hold off
grid on
ylabel("\deltae_y' [m]")

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
xlabel('Orbital Periods')

% figure
% subplot(3,2,1)
% hold on
% plot(0, ROE_PPM(1) * a_chief * 1e3, 'ro')
% plot(reconfiguration_tspan / T_chief, ROE_history(:,1) * a_chief * 1e3, 'b-')
% plot(reconfiguration_tspan(end) / T_chief, ROE_IAM(1) * a_chief * 1e3, 'rx')
% hold off
% legend('Start', 'Maneuvering', 'Target')
% grid on
% ylabel('\deltaa [m]')
% 
% subplot(3,2,3)
% hold on
% plot(0, ROE_PPM(2) * a_chief * 1e3, 'ro')
% plot(reconfiguration_tspan / T_chief, ROE_history(:,2) * a_chief * 1e3, 'b-')
% plot(reconfiguration_tspan(end) / T_chief, ROE_IAM(2) * a_chief * 1e3, 'rx')
% hold off
% grid on
% ylabel('\delta\lambda_e [m]')
% 
% subplot(3,2,5)
% hold on
% plot(0, ROE_PPM(3) * a_chief * 1e3, 'ro')
% plot(reconfiguration_tspan / T_chief, ROE_history(:,3) * a_chief * 1e3, 'b-')
% plot(reconfiguration_tspan(end) / T_chief, ROE_IAM(3) * a_chief * 1e3, 'rx')
% hold off
% grid on
% ylabel('\deltae_x [m]')
% xlabel('Orbital Periods')
% 
% subplot(3,2,2)
% hold on
% plot(0, ROE_PPM(4) * a_chief * 1e3, 'ro')
% plot(reconfiguration_tspan / T_chief, ROE_history(:,4) * a_chief * 1e3, 'b-')
% plot(reconfiguration_tspan(end) / T_chief, ROE_IAM(4) * a_chief * 1e3, 'rx')
% hold off
% grid on
% ylabel('\deltae_y [m]')
% 
% subplot(3,2,4)
% hold on
% plot(0, ROE_PPM(5) * a_chief * 1e3, 'ro')
% plot(reconfiguration_tspan / T_chief, ROE_history(:,5) * a_chief * 1e3, 'b-')
% plot(reconfiguration_tspan(end) / T_chief, ROE_IAM(5) * a_chief * 1e3, 'rx')
% hold off
% grid on
% ylabel('\deltai_x [m]')
% 
% subplot(3,2,6)
% hold on
% plot(0, ROE_PPM(6) * a_chief * 1e3, 'ro')
% plot(reconfiguration_tspan / T_chief, ROE_history(:,6) * a_chief * 1e3, 'b-')
% plot(reconfiguration_tspan(end) / T_chief, ROE_IAM(6) * a_chief * 1e3, 'rx')
% hold off
% grid on
% ylabel('\deltai_y [m]')
% xlabel('Orbital Periods')


figure
subplot(3,1,1)
hold on
plot(modified_ROE_PPM(2) * a_chief * 1e3, modified_ROE_PPM(1) * a_chief * 1e3, 'ro')
plot(modified_ROE_history(:,2) * a_chief * 1e3, modified_ROE_history(:,1) * a_chief * 1e3, 'b-')
plot(modified_ROE_IAM(2) * a_chief * 1e3, modified_ROE_IAM(1) * a_chief * 1e3, 'rx')
hold off
grid on
axis equal
legend('Start', 'Maneuvering', 'Target')
xlabel('\delta\lambda_e [m]')
ylabel('\deltaa [m]')

subplot(3,1,2)
hold on
plot(modified_ROE_PPM(3) * a_chief * 1e3, modified_ROE_PPM(4) * a_chief * 1e3, 'ro')
plot(modified_ROE_history(:,3) * a_chief * 1e3, modified_ROE_history(:,4) * a_chief * 1e3, 'b-')
plot(modified_ROE_IAM(3) * a_chief * 1e3, modified_ROE_IAM(4) * a_chief * 1e3, 'rx')
hold off
grid on
axis equal
xlabel("\deltae_x' [m]")
ylabel("\deltae_y' [m]")

subplot(3,1,3)
hold on
plot(modified_ROE_PPM(5) * a_chief * 1e3, modified_ROE_PPM(6) * a_chief * 1e3, 'ro')
plot(modified_ROE_history(:,5) * a_chief * 1e3, modified_ROE_history(:,6) * a_chief * 1e3, 'b-')
plot(modified_ROE_IAM(5) * a_chief * 1e3, modified_ROE_IAM(6) * a_chief * 1e3, 'rx')
hold off
grid on
axis equal
xlabel('\deltai_x [m]')
ylabel('\deltai_y [m]')


figure
subplot(3,2,1)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,1) * a_chief * 1e3)
ylabel('\Delta\deltaa [m]')

subplot(3,2,3)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,2) * a_chief * 1e3)
ylabel('\Delta\delta\lambda_e [m]')

subplot(3,2,5)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,3) * a_chief * 1e3)
ylabel("\Delta\deltae_x' [m]")
xlabel('Orbital Periods')

subplot(3,2,2)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,4) * a_chief * 1e3)
ylabel("\Delta\deltae_y' [m]")

subplot(3,2,4)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,5) * a_chief * 1e3)
ylabel('\Delta\deltai_x [m]')

subplot(3,2,6)
plot(reconfiguration_tspan / T_chief, control_tracking_error_history(:,6) * a_chief * 1e3)
ylabel('\Delta\deltai_y [m]')
xlabel('Orbital Periods')

figure
hold on
plot(reconfiguration_tspan / T_chief, total_dv(:,1), 'b-')
plot(reconfiguration_tspan / T_chief, total_dv(:,2), 'g-')
plot(reconfiguration_tspan / T_chief, sqrt(total_dv(:,1).^2 + total_dv(:,2).^2), 'r-')
hold off
grid on
legend('\DeltaV_T', '\DeltaV_N', 'Total \DeltaV')
ylabel('\DeltaV [m/s]')
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
    B = zeros(5,2);
    eta = sqrt(1-e^2);

    B(1,1) = 2 * (1+e*cos(nu)) / eta;

    B(2,1) = eta * (e + cos(nu)*(2+e*cos(nu))) / (1+e*cos(nu));

    B(3,1) = (eta / e) * sin(nu) * (2+e*cos(nu)) / (1+e*cos(nu));

    B(4,2) = eta * cos(omega + nu) / (1+e*cos(nu));

    B(5,2) = eta * sin(omega + nu) / (1+e*cos(nu));

    B = B ./ (a .* n);
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

function nu = inplane_planning(e, omega, i, delta_ex, delta_ey, delta_iy)
    syms x
    f = tan(x) * (1 + e * sin(omega) / (sin(x) * (2 + e * cos(omega) * cos(x) + e * sin(omega) * sin(x)))) / (1 + e * cos(omega) / (cos(x) * (2 + e * cos(omega) * cos(x) + e * sin(omega) * sin(x)))) - (delta_ey * tan(i) + e * cos(omega) * delta_iy) / (delta_ex * tan(i) - e * sin(omega) * delta_iy);
    nu = vpasolve(f,x,[0 2 * pi]);
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