%% Problem 2

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2

%% Part a

% For the chief
a_c = 36943; % km
e_c = 0.3;
inc_c = deg2rad(59);
omega_c = deg2rad(188);
RAAN_c = deg2rad(84);
nu_c = deg2rad(0);
M_c = 0;
T = 2 * pi * sqrt(a_c^3 / mu);
n_c = sqrt(mu / a_c^3);

% For the deputy
delta_a = 0;
delta_e = - e_c * 0.00001 * a_c;
delta_i = deg2rad(0.0005) * a_c;
delta_omega = 0;
delta_RAAN = - deg2rad(0.001) * a_c;
delta_nu = deg2rad(-0.00016) * a_c;

% Inertial position and velocity in ECI
a_d = a_c + delta_a;
e_d = e_c + delta_e / a_c;
inc_d = inc_c + delta_i / a_c;
omega_d = omega_c + delta_omega / a_c;
RAAN_d = RAAN_c + delta_RAAN / a_c;
nu_d = nu_c + delta_nu / a_c;
E_d = atan2(sin(nu_d)*sqrt(1-e_d^2)/(1+e_d*cos(nu_d)), (e_d+cos(nu_d))/(1+e_d*cos(nu_d)));
M_d = E_d - e_d * sin(E_d);

[pos_chief, vel_chief] = OE2ECI(a_c, e_c, inc_c, omega_c, RAAN_c, nu_c, mu);
[pos_deputy, vel_deputy] = OE2ECI(a_d, e_d, inc_d, omega_d, RAAN_d, nu_d, mu);

% Relative position and velocity in ECI
rho_ECI = pos_deputy - pos_chief;
rho_dot_ECI = vel_deputy - vel_chief;

% Relative position and velocity in RTN
state_RTN = ECI2RTN([pos_chief; vel_chief], [rho_ECI; rho_dot_ECI], mu);
rho_RTN = state_RTN(1:3)';
rho_dot_RTN = state_RTN(4:6)';

norm(rho_RTN) / (0.001 * a_c)

%% Part b:

k = 1 + e_c * cos(nu_c);
k_prime = - e_c * sin(nu_c);
eta = sqrt(1 - e_c^2);
scaling_matrix = [a_c * eta^2 * eye(3), zeros(3);
       zeros(3), a_c * n_c / eta * eye(3)];
int_const = (scaling_matrix * STM_YA(nu_c, e_c, 0, n_c))^(-1) * [rho_RTN; rho_dot_RTN];

M_2_inv = [2*k^2*(k+1)/eta^2, 2*k^2*k_prime/eta^2, 0, -2*k_prime/eta^2, 2*k/eta^2, 0;
           (1-(k+1)^2/eta^2)*sin(nu_c), -(k+1)*k_prime/eta^2*sin(nu_c), 0, 1/eta^2*(cos(nu_c)-2*e_c/k), -1/eta^2*(1+1/k)*sin(nu_c), 0;
           -k/eta^2*(2*e_c + (2+k)*cos(nu_c)), -k_prime/eta^2*(e_c + (k+1)*cos(nu_c)), 0, -sin(nu_c)/eta^2, -1/eta^2*(e_c/k + (1+1/k)*cos(nu_c)), 0;
           (k+1)^2*k_prime/eta^2, k/eta^2*(2+k-k^2)-1, 0, 1/eta^2*(k-1-2/k), k_prime/eta^2*(1+1/k), 0;
           0, 0, sin(nu_c), 0, 0, cos(nu_c)/k;
           0, 0, e_c+cos(nu_c), 0, 0, -sin(nu_c)/k];
K = M_2_inv * [1/ (a_c * eta^2) * eye(3), zeros(3); zeros(3), eta/(a_c*n_c)*eye(3)] * [rho_RTN; rho_dot_RTN];

% int_const(1) = 0;

% k*rho_dot_RTN(2) + e_c*sin(nu_c)*(rho_dot_RTN(1)-rho_RTN(2)) + (2+e_c*cos(nu_c))*rho_RTN(1)
% % rho_RTN(1) = - (k*rho_dot_RTN(2) + e_c*sin(nu_c)*(rho_dot_RTN(1)-rho_RTN(2))) / (2+e_c*cos(nu_c));
% % rho_dot_RTN(2) = - (e_c*sin(nu_c)*(rho_dot_RTN(1)-rho_RTN(2)) + (2+e_c*cos(nu_c))*rho_RTN(1)) / k;
% k*rho_dot_RTN(2) + e_c*sin(nu_c)*(rho_dot_RTN(1)-rho_RTN(2)) + (2+e_c*cos(nu_c))*rho_RTN(1)

%% Part c: propagation of YA solution

N = 10000;
tspan = linspace(0, 15 * T, N);
relative_motion = zeros(N, 6);
relative_motion_2 = zeros(N, 6);

for i=1:N
    M = wrapTo2Pi(n_c * tspan(i));
    E = M2E(M, e_c, 1e-12);
    nu = atan2(sin(E)*sqrt(1-e_c^2)/(1-e_c*cos(E)), (cos(E)-e_c)/(1-e_c*cos(E)));
    STM = STM_YA(nu, e_c, tspan(i), n_c);
    relative_motion(i, :) = scaling_matrix * STM * int_const;
    relative_motion_2(i, :) = scaling_matrix * STM * K;
end


% Maybe put a star at (0,0,0) to indicate where the chief is?
% figure
% subplot(2,2,1)
% plot(relative_motion(:, 2), relative_motion(:, 1))
% grid on
% axis equal
% xlabel('T-axis')
% ylabel('R-axis')
% 
% subplot(2,2,2)
% plot(relative_motion(:, 3), relative_motion(:, 1))
% grid on
% axis equal
% xlabel('N-axis')
% ylabel('R-axis')
% 
% subplot(2,2,3)
% plot(relative_motion(:, 2), relative_motion(:, 3))
% grid on
% axis equal
% xlabel('T-axis')
% ylabel('N-axis')
% 
% subplot(2,2,4)
% plot3(relative_motion(:, 1), relative_motion(:, 2), relative_motion(:, 3))
% grid on
% axis equal
% view(3)
% xlabel('R-axis')
% ylabel('T-axis')
% zlabel('N-axis')

figure
subplot(2,2,1)
plot(relative_motion_2(:, 2), relative_motion_2(:, 1))
grid on
axis equal
xlabel('T-axis')
ylabel('R-axis')

subplot(2,2,2)
plot(relative_motion_2(:, 3), relative_motion_2(:, 1))
grid on
axis equal
xlabel('N-axis')
ylabel('R-axis')

subplot(2,2,3)
plot(relative_motion_2(:, 2), relative_motion_2(:, 3))
grid on
axis equal
xlabel('T-axis')
ylabel('N-axis')

subplot(2,2,4)
plot3(relative_motion_2(:, 1), relative_motion_2(:, 2), relative_motion_2(:, 3))
grid on
axis equal
view(3)
xlabel('R-axis')
ylabel('T-axis')
zlabel('N-axis')

figure
subplot(2,2,1)
plot(relative_motion_2(:, 5) / a_c, relative_motion_2(:, 4) / a_c)
grid on
axis equal
xlabel('T-axis')
ylabel('R-axis')

subplot(2,2,2)
plot(relative_motion_2(:, 6) / a_c, relative_motion_2(:, 4) / a_c)
grid on
axis equal
xlabel('N-axis')
ylabel('R-axis')

subplot(2,2,3)
plot(relative_motion_2(:, 5) / a_c, relative_motion_2(:, 6) / a_c)
grid on
axis equal
xlabel('T-axis')
ylabel('N-axis')

subplot(2,2,4)
plot3(relative_motion_2(:, 4) / a_c, relative_motion_2(:, 5) / a_c, relative_motion_2(:, 6) / a_c)
grid on
axis equal
view(3)
xlabel('R-axis')
ylabel('T-axis')
zlabel('N-axis')

%% Part e: Relative Orbital Elements

ROE = [(a_d - a_c) / a_c;
        M_d + omega_d - M_c - omega_c + (RAAN_d - RAAN_c) * cos(inc_c);
        e_d * cos(omega_d) - e_c * cos(omega_c);
        e_d * sin(omega_d) - e_c * sin(omega_c);
        inc_d - inc_c;
        (RAAN_d - RAAN_c) * sin(inc_c)]

%% Part f: Propagation using linear mapping

relative_motion_ROE = zeros(N,6);

for i=1:N
    M = wrapTo2Pi(n_c * tspan(i));
    E = M2E(M, e_c, 1e-15);
    nu = atan2(sin(E)*sqrt(1-e_c^2)/(1-e_c*cos(E)), (cos(E)-e_c)/(1-e_c*cos(E)));
    STM = STM_ROE(nu, e_c, omega_c, inc_c, tspan(i), n_c);
    relative_motion_ROE(i, :) = scaling_matrix * STM * ROE;
end

figure
subplot(2,2,1)
plot(relative_motion_ROE(:, 2), relative_motion_ROE(:, 1))
grid on
axis equal
xlabel('T-axis')
ylabel('R-axis')

subplot(2,2,2)
plot(relative_motion_ROE(:, 3), relative_motion_ROE(:, 1))
grid on
axis equal
xlabel('N-axis')
ylabel('R-axis')

subplot(2,2,3)
plot(relative_motion_ROE(:, 2), relative_motion_ROE(:, 3))
grid on
axis equal
xlabel('T-axis')
ylabel('N-axis')

subplot(2,2,4)
plot3(relative_motion_ROE(:, 1), relative_motion_ROE(:, 2), relative_motion_ROE(:, 3))
grid on
axis equal
view(3)
xlabel('R-axis')
ylabel('T-axis')
zlabel('N-axis')

figure
subplot(2,2,1)
plot(relative_motion_ROE(:, 5) / a_c, relative_motion_ROE(:, 4) / a_c)
grid on
axis equal
xlabel('T-axis')
ylabel('R-axis')

subplot(2,2,2)
plot(relative_motion_ROE(:, 6) / a_c, relative_motion_ROE(:, 4) / a_c)
grid on
axis equal
xlabel('N-axis')
ylabel('R-axis')

subplot(2,2,3)
plot(relative_motion_ROE(:, 5) / a_c, relative_motion_ROE(:, 6) / a_c)
grid on
axis equal
xlabel('T-axis')
ylabel('N-axis')

subplot(2,2,4)
plot3(relative_motion_ROE(:, 4) / a_c, relative_motion_ROE(:, 5) / a_c, relative_motion_ROE(:, 6) / a_c)
grid on
axis equal
view(3)
xlabel('R-axis')
ylabel('T-axis')
zlabel('N-axis')

%% Part g: Relationship between YA integration constants and ROEs

rot = [1, 0, 0, 0, 0, 0;
       0, cos(omega_c), sin(omega_c), 0, 0, 0;
       0, -sin(omega_c), cos(omega_c), 0, 0, 0;
       0, 0, 0, 1, 0, 0;
       0, 0, 0, 0, cos(omega_c), sin(omega_c);
       0, 0, 0, 0, -sin(omega_c), cos(omega_c)];
K_u = rot * K;

e_x = e_c * cos(omega_c);
e_y = e_c * sin(omega_c);

mat = [1, 0, 0, 0, 0, 0;
       0, -e_x*(eta+1/(1+eta)), e_y*(eta+1/(1+eta)), 1, 0, 0;
       0, e_x*e_y, e_x^2-1, -e_y, 0, -e_y*cot(inc_c);
       0, e_y^2-1, e_x*e_y, e_x, 0, e_x*cot(inc_c);
       0, 0, 0, 0, 1, 0;
       0, 0, 0, 0, 0, -1];
ROE_from_YA = mat * K_u;
ROE_from_YA - ROE

%% Part h: Comparison with numerical integration of FODE

IC = [pos_chief', vel_chief', pos_deputy', vel_deputy'];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100);
[t_FODE, y_FODE] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, IC, options);

relative_motion_FODE = zeros(N, 6);
errors_YA = zeros(N, 6);
errors_ROE = zeros(N, 6);

for i=1:N
    rel_state_ECI = y_FODE(i, 7:12) - y_FODE(i, 1:6);
    rel_state_RTN = ECI2RTN(y_FODE(i, 1:6), rel_state_ECI, mu);
    relative_motion_FODE(i, :) = rel_state_RTN;
%     abs(rel_state_RTN - relative_motion_2(i, :))
%     abs(rel_state_RTN - relative_motion_ROE(i, :))
    errors_YA(i, :) = abs(rel_state_RTN - relative_motion_2(i, :));
    errors_ROE(i, :) = abs(rel_state_RTN - relative_motion_ROE(i, :));
end

figure
subplot(2,2,1)
plot(relative_motion_FODE(:, 2), relative_motion_FODE(:, 1))
grid on
axis equal
xlabel('T-axis')
ylabel('R-axis')

subplot(2,2,2)
plot(relative_motion_FODE(:, 3), relative_motion_FODE(:, 1))
grid on
axis equal
xlabel('N-axis')
ylabel('R-axis')

subplot(2,2,3)
plot(relative_motion_FODE(:, 2), relative_motion_FODE(:, 3))
grid on
axis equal
xlabel('T-axis')
ylabel('N-axis')

subplot(2,2,4)
plot3(relative_motion_FODE(:, 1), relative_motion_FODE(:, 2), relative_motion_FODE(:, 3))
grid on
axis equal
view(3)
xlabel('R-axis')
ylabel('T-axis')
zlabel('N-axis')

figure
subplot(2,2,1)
plot(relative_motion_FODE(:, 5) / a_c, relative_motion_FODE(:, 4) / a_c)
grid on
axis equal
xlabel('T-axis')
ylabel('R-axis')

subplot(2,2,2)
plot(relative_motion_FODE(:, 6) / a_c, relative_motion_FODE(:, 4) / a_c)
grid on
axis equal
xlabel('N-axis')
ylabel('R-axis')

subplot(2,2,3)
plot(relative_motion_FODE(:, 5) / a_c, relative_motion_FODE(:, 6) / a_c)
grid on
axis equal
xlabel('T-axis')
ylabel('N-axis')

subplot(2,2,4)
plot3(relative_motion_FODE(:, 4) / a_c, relative_motion_FODE(:, 5) / a_c, relative_motion_FODE(:, 6) / a_c)
grid on
axis equal
view(3)
xlabel('R-axis')
ylabel('T-axis')
zlabel('N-axis')

figure
subplot(3,2,1)
title('Error in position')
hold on
plot(tspan / T, errors_YA(:, 1))
plot(tspan / T, errors_ROE(:, 1))
hold off
legend('YA','ROE')
grid on
ylabel('R-direction [km]')

subplot(3,2,3)
hold on
plot(tspan / T, errors_YA(:, 2))
plot(tspan / T, errors_ROE(:, 2))
hold off
grid on
ylabel('T-direction [km]')

subplot(3,2,5)
hold on
plot(tspan / T, errors_YA(:, 3))
plot(tspan / T, errors_ROE(:, 3))
hold off
grid on
ylabel('N-direction [km]')
xlabel('Orbital Periods')

subplot(3,2,2)
title('Error in velocity')
hold on
plot(tspan / T, errors_YA(:, 4))
plot(tspan / T, errors_ROE(:, 4))
hold off
grid on
ylabel('R-direction [km/s]')

subplot(3,2,4)
hold on
plot(tspan / T, errors_YA(:, 5))
plot(tspan / T, errors_ROE(:, 5))
hold off
grid on
ylabel('T-direction [km/s]')

subplot(3,2,6)
hold on
plot(tspan / T, errors_YA(:, 6))
plot(tspan / T, errors_ROE(:, 6))
hold off
grid on
ylabel('N-direction [km/s]')


%% Functions

function M = STM_YA(f, e, t, n)
    k = 1 + e * cos(f);
    k_prime = - e * sin(f);
    eta = sqrt(1 - e^2);
    tau = n * t / eta^3;
    
    M = [      1/k + 3/2*k_prime*tau,          sin(f),            cos(f),       0,          0,          0;
                            -3/2*k*tau,  (1+1/k)*cos(f),   -(1+1/k)*sin(f),     1/k,          0,          0;
                                     0,               0,                 0,       0, 1/k*sin(f), 1/k*cos(f);
           k_prime/2-3/2*k^2*(k-1)*tau,      k^2*cos(f),       -k^2*sin(f),       0,          0,          0;
               -3/2*(k+k^2*k_prime*tau), -(k^2+1)*sin(f), -e-(k^2+1)*cos(f), -k_prime,          0,          0;
                                     0,               0,                 0,       0,   e+cos(f),    -sin(f)];
end

function M = STM_ROE(f, e, omega, i, t, n)
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