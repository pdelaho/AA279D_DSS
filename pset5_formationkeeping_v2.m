%% For station keeping

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
R_E = 6378; % km
J2 = 0.108263e-2;
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100); % options for numerical integration

% Be careful between mean and osculating elements this time, don't forget
% to set tolerance = 1e-12 when going from mean to true anomaly

% First just integrate the ROE over the 6 hours of formation keeping to see
% what the boundaries of the keep in zone should be to have a few maneuvers
% to do

% Defining the ROE for Inertial Attitude Mode (IAM)

% Initial conditions in mean elements for the chief

a_chief = 36943; % km
e_chief = 0.8111;
inc_chief = deg2rad(59);
omega_chief = deg2rad(188);
RAAN_chief = deg2rad(84);
n_chief = sqrt(mu / a_chief^3);
T_chief = 2 * pi / n_chief;
M_chief = n_chief * (6 * 3600 + 49 * 60);
% nu_chief_IAM = mean2true(M_chief, e_chief);
oe_mean_chief_IAM = [a_chief, e_chief, inc_chief, omega_chief, RAAN_chief, M_chief];
% temp = mean2osc([oe_mean_chief_IAM(1)*1e3 oe_mean_chief_IAM(2:3) oe_mean_chief_IAM(5) oe_mean_chief_IAM(4) oe_mean_chief_IAM(6)]);
% oe_osc_chief_IAM = [temp(1)*1e-3 temp(2:3)' temp(5) temp(4) temp(6)];
% nu_osc_chief_IAM = mean2true(temp(6), temp(2));

% Initial conditions in mean elements for the deputy

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
nu_deputy_IAM = mean2true(M_deputy, e_deputy, 1e-12);
oe_mean_deputy_IAM = [a_deputy, e_deputy, inc_deputy, omega_deputy, RAAN_deputy, M_deputy];
% temp = mean2osc([oe_mean_deputy_IAM(1)*1e3 oe_mean_deputy_IAM(2:3) oe_mean_deputy_IAM(5) oe_mean_deputy_IAM(4) oe_mean_deputy_IAM(6)]);
% oe_osc_deputy_IAM = [temp(1)*1e-3 temp(2:3)' temp(5) temp(4) temp(6)];
% nu_osc_deputy_IAM = mean2true(temp(6), temp(2), 1e-12);
[pos_deputy_IAM, vel_deputy_IAM] = OE2ECI(oe_mean_deputy_IAM(1), oe_mean_deputy_IAM(2), oe_mean_deputy_IAM(3), oe_mean_deputy_IAM(4), oe_mean_deputy_IAM(5), nu_deputy_IAM, mu);

ROE_IAM = [delta_a, delta_lambda, delta_ex, delta_ey, delta_ix, delta_iy]';

init_cond = [oe_mean_chief_IAM, oe_mean_deputy_IAM];
N = 10000;
tspan = linspace(6 * 3600 + 49 * 20, 12 * 3600 + 49 * 60, N);
% tspan = linspace(0, T_chief, N);
dt = tspan(2) - tspan(1);
[t, y] = ode89(@(t, state) GVE(t, state, mu, J2, R_E), tspan, init_cond, options);

ROE_history = zeros(N, 6);
rel_dist = zeros(N,1);
for j=1:N
    oe_chief = y(j, 1:6); % a, e, i, omega, RAAN, M
    oe_deputy = y(j, 7:12);

    ROE_history(j, :) = OE2ROE(oe_chief, oe_deputy);
    [pos_deputy, vel_deputy] = OE2ECI(oe_deputy(1), oe_deputy(2), oe_deputy(3), oe_deputy(4), oe_deputy(5), mean2true(oe_deputy(6), oe_deputy(2), 1e-14), mu);
    rel_dist(j) = norm(pos_deputy - pos_deputy_IAM);
end

% trying to do station keeping

c = [0, 80e-3 / a_chief, -20e-3 / a_chief, 89.4e-3 / a_chief, 0, 79e-3 / a_chief];
R = [5e-6 / a_chief, 5e-6 / a_chief, 1.5e-5 / a_chief];
eps = 1e-26;

prev_oe_chief = oe_mean_chief_IAM;
prev_oe_deputy = oe_mean_deputy_IAM;
prev_t = 6 * 3600 + 49 * 60;
ROE = [ROE_IAM'];
time = [6 * 3600 + 49 * 60];
j = 1;

while j <= N && prev_t < tspan(end)
    j
%     j
    % check if the deputy is too far from the desired state
    % assuming we have prev_oe_chief and prev_oe_deputy in terms of
    % osculating elements

    prev_ROE = OE2ROE(prev_oe_chief, prev_oe_deputy); % should be in terms of mean elements
    (prev_ROE - ROE_IAM') * a_chief * 1e3
    (prev_ROE(1) - c(1))^2 + (prev_ROE(2) - c(2))^2 
    R(1)^2
%     prev_t

    if (prev_ROE(1) - c(1))^2 + (prev_ROE(2) - c(2))^2 > R(1)^2 + eps || (prev_ROE(1) - c(1))^2 + (prev_ROE(2) - c(2))^2 < R(1)^2 - eps || (prev_ROE(3) - c(3))^2 + (prev_ROE(4) - c(4))^2 > R(2)^2 + eps || (prev_ROE(3) - c(3))^2 + (prev_ROE(4) - c(4))^2 < R(2)^2 - eps || (prev_ROE(5) - c(5))^2 + (prev_ROE(6) - c(6))^2 > R(3)^2 + eps || (prev_ROE(5) - c(5))^2 + (prev_ROE(6) - c(6))^2 < R(3)^2 - eps
        % outside the keep-in-zone, need to maneuver, 6 maneuvers, one
        % every 5 minutes?
        % don't get out of this loop before the maneuvers are all executed
        % Apply the delta-v to the mean elements, not the osculating ones
        nb_steps = round(5 * 60 / (dt)); % how many time steps between maneuvers so that they are about 5 minutes apart
%         A_f0 = ROE_STM(n_chief, 6 * nb_steps * dt);
        A_f0 = STM_ROE_J2(6 * nb_steps * dt, prev_oe_chief(1), prev_oe_chief(2), prev_oe_chief(3), prev_oe_chief(4), J2, R_E, mu);
        B_0 = control_input(prev_oe_chief(1), sqrt(mu / prev_oe_chief(1)^3), prev_oe_chief(2), mean2true(prev_oe_chief(6), prev_oe_chief(2), 1e-12), prev_oe_chief(4), prev_oe_chief(3)); % a, n, e, nu, omega, i
        delta_ROE = ROE_IAM - A_f0 * prev_ROE';
        
%         A_f1 = ROE_STM(n_chief, 5 * nb_steps * dt);
        % should integrate the chief's mean elements between the different
        % time steps to be accurate
        span = linspace(prev_t, prev_t + nb_steps * dt, 2);
        [t, y] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), span, [prev_oe_chief, prev_oe_deputy], options);
        B_1 = control_input(y(end, 1), sqrt(mu / y(end, 1)^3), y(end, 2), mean2true(y(end, 6), y(end, 2), 1e-12), y(end, 4), y(end, 3));
        A_f1 = STM_ROE_J2(5 * nb_steps * dt, y(end,1), y(end,2), y(end,3), y(end,4), J2, R_E, mu);

%         A_f2 = ROE_STM(n_chief, 4 * nb_steps * dt);
        span = linspace(prev_t, prev_t + 2 * nb_steps * dt, 2);
        [t, y] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), span, [prev_oe_chief, prev_oe_deputy], options);
        B_2 = control_input(y(end, 1), sqrt(mu / y(end, 1)^3), y(end, 2), mean2true(y(end, 6), y(end, 2), 1e-12), y(end, 4), y(end, 3));
        A_f2 = STM_ROE_J2(4 * nb_steps * dt, y(end,1), y(end,2), y(end,3), y(end,4), J2, R_E, mu);

%         A_f3 = ROE_STM(n_chief, 3 * nb_steps * dt);
        span = linspace(prev_t, prev_t + 3 * nb_steps * dt, 2);
        [t, y] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), span, [prev_oe_chief, prev_oe_deputy], options);
        B_3 = control_input(y(end, 1), sqrt(mu / y(end, 1)^3), y(end, 2), mean2true(y(end, 6), y(end, 2), 1e-12), y(end, 4), y(end, 3));
        A_f3 = STM_ROE_J2(3 * nb_steps * dt, y(end,1), y(end,2), y(end,3), y(end,4), J2, R_E, mu);

%         A_f4 = ROE_STM(n_chief, 2 * nb_steps * dt);
        span = linspace(prev_t, prev_t + 4 * nb_steps * dt, 2);
        [t, y] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), span, [prev_oe_chief, prev_oe_deputy], options);
        B_4 = control_input(y(end, 1), sqrt(mu / y(end, 1)^3), y(end, 2), mean2true(y(end, 6), y(end, 2), 1e-12), y(end, 4), y(end, 3));
        A_f4 = STM_ROE_J2(2 * nb_steps * dt, y(end,1), y(end,2), y(end,3), y(end,4), J2, R_E, mu);

%         A_f5 = ROE_STM(n_chief, nb_steps * dt);
        span = linspace(prev_t, prev_t + 5 * nb_steps * dt, 2);
        [t, y] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), span, [prev_oe_chief, prev_oe_deputy], options);
        B_5 = control_input(y(end, 1), sqrt(mu / y(end, 1)^3), y(end, 2), mean2true(y(end, 6), y(end, 2), 1e-12), y(end, 4), y(end, 3));
        A_f5 = STM_ROE_J2(nb_steps * dt, y(end,1), y(end,2), y(end,3), y(end,4), J2, R_E, mu);

        G = [A_f0 * B_0, A_f1 * B_1, A_f2 * B_2, A_f3 * B_3, A_f4 * B_4, A_f5 * B_5];
        delta_v = pinv(G) * delta_ROE

        % for each delta-v, transform it from RTN to OE then integrate til
        % the next maneuver

        delta_ROE_0 = transform(prev_oe_chief, delta_v(1:3), mu);
        ic_01 = [prev_oe_chief, prev_oe_deputy + delta_ROE_0];
        tspan_01 = linspace(prev_t, prev_t + nb_steps * dt, 10);
        [t_01, y_01] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), tspan_01, ic_01, options);
        
%         ROE_01 = zeros(1000, 6);
        for k=1:10
            ROE = [ROE; OE2ROE(y_01(k, 1:6), y_01(k, 7:12))];
            time = [time, t_01(k)];
        end
        
        delta_ROE_1 = transform(y_01(end, 1:6), delta_v(4:6), mu);
        ic_12 = [y_01(end, 1:6), y_01(end, 7:12) + delta_ROE_1];
        tspan_12 = linspace(prev_t + nb_steps * dt, prev_t + 2 * nb_steps * dt, 10);
        [t_12, y_12] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), tspan_12, ic_12, options);
        
%         ROE_12 = zeros(1000, 6);
        for k=1:10
            ROE = [ROE; OE2ROE(y_12(k, 1:6), y_12(k, 7:12))];
            time = [time, t_12(k)];
        end
        
        delta_ROE_2 = transform(y_12(end, 1:6), delta_v(7:9), mu);
        ic_23 = [y_12(end, 1:6), y_12(end, 7:12) + delta_ROE_2];
        tspan_23 = linspace(prev_t + 2 * nb_steps * dt, prev_t + 3 * nb_steps * dt, 10);
        [t_23, y_23] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), tspan_23, ic_23, options);
        
%         ROE_23 = zeros(1000, 6);
        for k=1:10
            ROE = [ROE; OE2ROE(y_23(k, 1:6), y_23(k, 7:12))];
            time = [time, t_23(k)];
        end
        
        delta_ROE_3 = transform(y_23(end, 1:6), delta_v(10:12), mu);
        ic_34 = [y_23(end, 1:6), y_23(end, 7:12) + delta_ROE_3];
        tspan_34 = linspace(prev_t + 3 * nb_steps * dt, prev_t + 4 * nb_steps * dt, 10);
        [t_34, y_34] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), tspan_34, ic_34, options);
        
%         ROE_34 = zeros(1000, 6);
        for k=1:10
            ROE = [ROE; OE2ROE(y_34(k, 1:6), y_34(k, 7:12))];
            time = [time, t_34(k)];
        end
        
        delta_ROE_4 = transform(y_34(end, 1:6), delta_v(13:15), mu);
        ic_45 = [y_34(end, 1:6), y_34(end, 7:12) + delta_ROE_4];
        tspan_45 = linspace(prev_t + 4 * nb_steps * dt, prev_t + 5 * nb_steps * dt, 10);
        [t_45, y_45] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), tspan_45, ic_45, options);
        
%         ROE_45 = zeros(1000, 6);
        for k=1:10
            ROE = [ROE; OE2ROE(y_45(k, 1:6), y_45(k, 7:12))];
            time = [time, t_45(k)];
        end
        
        delta_ROE_5 = transform(y_45(end, 1:6), delta_v(16:18), mu);
        ic_56 = [y_45(end, 1:6), y_45(end, 7:12) + delta_ROE_5];
        tspan_56 = linspace(prev_t + 5 * nb_steps * dt, prev_t + 6 * nb_steps * dt, 10);
        [t_56, y_56] = ode89(@ (t, state) GVE(t, state, mu, J2, R_E), tspan_56, ic_56, options);
        
%         ROE_56 = zeros(1000, 6);
        for k=1:10
            ROE = [ROE; OE2ROE(y_56(k, 1:6), y_56(k, 7:12))];
            time = [time, t_56(k)];
        end

        prev_oe_chief = y_56(end, 1:6);
        prev_oe_deputy = y_56(end, 7:12);
        prev_t = t_56(end);
        % end of the maneuvers
        j = j + 1;
    else
        % no need for maneuvers, just integrate for the next time step 
        span = linspace(prev_t, prev_t + dt, 10);
        ic = [prev_oe_chief, prev_oe_deputy];
        [t, y] = ode89(@(t, state) GVE(t, state, mu, J2, R_E), span, ic, options);

        for k=1:10
            ROE = [ROE; OE2ROE(y(k, 1:6), y(k,7:12))];
            time = [time, t(k)];
            prev_oe_chief = y(end, 1:6);
            prev_oe_deputy = y(end, 7:12);
            prev_t = t(end);
        end
        j = j + 1;
    end
    
end

%% Plotting

figure
subplot(3,2,1)
plot((ROE_history(:, 1) - ROE_IAM(1)) * a_chief * 1e3)
grid on
ylabel('\deltaa [m]')

subplot(3,2,3)
plot((ROE_history(:, 2) - ROE_IAM(2)) * a_chief * 1e3)
grid on
ylabel('\delta\lambda [m]')

subplot(3,2,5)
plot((ROE_history(:, 3) - ROE_IAM(3)) * a_chief * 1e3)
grid on
ylabel('\deltae_x [m]')

subplot(3,2,2)
plot((ROE_history(:, 4) - ROE_IAM(4)) * a_chief * 1e3)
grid on
ylabel('\deltae_y [m]')

subplot(3,2,4)
plot((ROE_history(:, 5) - ROE_IAM(5)) * a_chief * 1e3)
grid on
ylabel('\deltai_x [m]')

subplot(3,2,6)
plot((ROE_history(:, 6) - ROE_IAM(6)) * a_chief * 1e3)
grid on
ylabel('\deltai_y [m]')

figure
plot(tspan / T_chief, rel_dist)


%% Functions

function ROE = OE2ROE(oe_chief, oe_deputy)
    ROE = zeros(1, 6);
    ROE(1) = (oe_deputy(1) - oe_chief(1)) / oe_chief(1);
    ROE(2) = wrapToPi((oe_deputy(6) + oe_deputy(4)) - (oe_chief(6) + oe_chief(4)) + (oe_deputy(5) - oe_chief(5)) * cos(oe_chief(3)));
    ROE(3) = oe_deputy(2) * cos(oe_deputy(4)) - oe_chief(2) * cos(oe_chief(4));
    ROE(4) = oe_deputy(2) * sin(oe_deputy(4)) - oe_chief(2) * sin(oe_chief(4));
    ROE(5) = oe_deputy(3) - oe_chief(3);
    ROE(6) = wrapToPi((oe_deputy(5) - oe_chief(5)) * sin(oe_chief(3)));
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

function A = ROE_STM(t, n) % not the right STM, need to use the one with the J2 effects
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

function statedot = GVE(t, state, mu, J2, R_E)
    % for 2 satellites
    statedot = zeros(size(state));
    n_c = sqrt(mu / state(1)^3);
    statedot(4) = 3 / 4 * n_c * J2 * (R_E / (state(1) * (1 - state(2)^2)))^2 * (5 * cos(state(3))^2 - 1);
    statedot(5) = -3 / 2 * n_c * J2 * (R_E / (state(1) * (1 - state(2)^2)))^2 * cos(state(3));
    statedot(6) = n_c + 3 / 4 * n_c * J2 * (R_E / (state(1) * (1 - state(2)^2)))^2 * sqrt(1 - state(2)^2) * (3 * cos(state(3))^2 - 1);
    n_d = sqrt(mu / state(7)^3);
    statedot(10) = 3 / 4 * n_d * J2 * (R_E / (state(7) * (1 - state(8)^2)))^2 * (5 * cos(state(9))^2 - 1);
    statedot(11) = -3 / 2 * n_d * J2 * (R_E / (state(7) * (1 - state(8)^2)))^2 * cos(state(9));
    statedot(12) = n_d + 3 / 4 * n_d * J2 * (R_E / (state(7) * (1 - state(8)^2)))^2 * sqrt(1 - state(8)^2) * (3 * cos(state(9))^2 - 1);
    % a, e, i, omega, RAAN, M 
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

