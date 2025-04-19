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

[pos_chief, vel_chief] = OEtoECI(a_c, e_c, inc_c, omega_c, RAAN_c, nu_c, mu);
[pos_deputy, vel_deputy] = OEtoECI(a_d, e_d, inc_d, omega_d, RAAN_d, nu_d, mu);

% Relative position and velocity in RTN
rho_ECI = pos_deputy - pos_chief;
rho_dot_ECI = vel_deputy - vel_chief;

rot_ECItoRTN = ECItoRTN([pos_chief; vel_chief]);
omega_vec = [0 0 sqrt(mu / (a_c^3 * (1 - e_c^2)^3)) * (1 + e_c * cos(nu_c))^2];

rho_RTN = rot_ECItoRTN * rho_ECI;
rho_dot_RTN = rot_ECItoRTN * rho_dot_ECI - cross(omega_vec', rho_RTN);

norm(rho_RTN) / (0.001 * a_c)

%% Part b:

% Do the deput and chief need to have the same eccentricity?
% If not, which e is it in the equations? The chief?

% k = 1 + e_c * cos(nu_c);
% x1_bar = k * sin(nu_c);
% cosE = (e_c + cos(nu_c)) / k;
% sinE = sin(nu_c) * sqrt(1 - e_c^2) / k;
% E = atan2(sinE,cosE);
% % what values is E supposed to be between?
% 
% K1 = (1 - e_c^2)^(-5/2) * (- 3 * e_c * E / 2 + (1 + e_c^2) * sinE - e_c / 4 * sin(2*E));
% % x2_bar = 2 * e_c * x1_bar * K1 - cos(nu_c) / k;
% % x3_bar = -2 * k * sin(nu_c) * K1 - cos(nu_c)^2 / k - cos(nu_c)^2;
% K2 = (1 - e_c^2)^(-5/2) * (1 / 2 * E - 1 / 4 * sin(2 * E) - e_c / 3 * sinE^3);
% x2_bar = 2 * e_c * x1_bar * (sin(nu_c) / k^3 - 3 * e_c * K2) - cos(nu_c) / k;
% % Which formula for x2_bar and x3_bar to use?
% x3_bar = 6 * e_c * x1_bar * K2 - 2 * sin(nu_c)^2 / k^2 - cos(nu_c)^2 / k - cos(nu_c)^2;
% S1 = -cos(nu_c) * (1 + e_c / 2 * cos(nu_c));
% S2 = 3 * e_c * k^2 * K2 - sin(nu_c) / k;
% S3 = - 6 * k^2 * K2 - (2 + k) / (2 * k) * sin(2 * nu_c);

% c = k * cos(nu_c);
% s = k * sin(nu_c);
% phi_inv = 1 / eta^2 * [-3*s*(k+e_c^2)/k^2, c-2*e_c, 0, -s*(k+1)/k, 0, 0;
%                        -3*(e_c+c/k), -s, 0, -(c*(k+1)/k)+e_c, 0, 0;
%                        3*k-eta^2, e_c*s, 0, k^2, 0, 0;
%                        -3*e_c*s*(k+1)/k^2, -2+e_c*c, eta^2, -e_c*s*(k+1)/k, 0, 0;
%                        0, 0, 0, 0, eta^2*cos(nu_c), -eta^2*sin(nu_c);
%                        0, 0, 0, 0, eta^2*sin(nu_c), eta^2*cos(nu_c)];
% % Computing integration constants
% int_const = phi_inv * [rho_RTN; rho_dot_RTN]

k = 1 + e_c * cos(nu_c);
k_prime = - e_c * sin(nu_c);
eta = sqrt(1 - e_c^2);
M_1 = [a_c * eta^2 * eye(3), zeros(3);
       zeros(3), a_c * n_c / eta * eye(3)]
int_const = (STM_YA(nu_c, e_c, 0, n_c, M_1))^(-1) * [rho_RTN; rho_dot_RTN]

M_2_inv = [2*k^2*(k+1)/eta^2, 2*k^2*k_prime/eta^2, 0, -2*k_prime/eta^2, 2*k/eta^2, 0;
           (1-(k+1)^2/eta^2)*sin(nu_c), -(k+1)*k_prime/eta^2*sin(nu_c), 0, 1/eta^2*(cos(nu_c)-2*e_c/k), -1/eta^2*(1+1/k)*sin(nu_c), 0;
           -k/eta^2*(2*e_c + (2+k)*cos(nu_c)), -k_prime/eta^2*(e_c + (k+1)*cos(nu_c)), 0, -sin(nu_c)/eta^2, -1/eta^2*(e_c/k + (1+1/k)*cos(nu_c)), 0;
           (k+1)^2*k_prime/eta^2, k/eta^2*(2+k-k^2)-1, 0, 1/eta^2*(k-1-2/k), k_prime/eta^2*(1+1/k), 0;
           0, 0, sin(nu_c), 0, 0, cos(nu_c)/k;
           0, 0, e_c+cos(nu_c), 0, 0, -sin(nu_c)/k];
int_const_2 = M_2_inv * [1/ (a_c * eta^2) * eye(3), zeros(3); zeros(3), eta/(a_c*n_c)*eye(3)] * [rho_RTN; rho_dot_RTN]

% int_const(1) = 0;

k*rho_dot_RTN(2) + e_c*sin(nu_c)*(rho_dot_RTN(1)-rho_RTN(2)) + (2+e_c*cos(nu_c))*rho_RTN(1)
% rho_RTN(1) = - (k*rho_dot_RTN(2) + e_c*sin(nu_c)*(rho_dot_RTN(1)-rho_RTN(2))) / (2+e_c*cos(nu_c));
% rho_dot_RTN(2) = - (e_c*sin(nu_c)*(rho_dot_RTN(1)-rho_RTN(2)) + (2+e_c*cos(nu_c))*rho_RTN(1)) / k;
k*rho_dot_RTN(2) + e_c*sin(nu_c)*(rho_dot_RTN(1)-rho_RTN(2)) + (2+e_c*cos(nu_c))*rho_RTN(1)

%% Part c: propagation of YA solution

N = 10000;
tspan = linspace(0, 15 * T, N);
relative_motion = zeros(N, 6);
relative_motion_2 = zeros(N, 6);

for i=1:N
    M = wrapTo2Pi(n_c * tspan(i));
    E = M2E(M, e_c, 1e-12);
    nu = atan2(sin(E)*sqrt(1-e_c^2)/(1-e_c*cos(E)), (cos(E)-e_c)/(1-e_c*cos(E)));
    STM = STM_YA(nu, e_c, tspan(i), n_c, M_1);
    relative_motion(i, :) = STM * int_const;
    relative_motion_2(i, :) = STM * int_const_2;
end

STM_YA(pi/4,e_c,10,n_c,M_1)


% Maybe put a star at (0,0,0) to indicate where the chief is?
figure
subplot(2,2,1)
plot(relative_motion(:, 2), relative_motion(:, 1))
grid on
axis equal
xlabel('T-axis')
ylabel('R-axis')

subplot(2,2,2)
plot(relative_motion(:, 3), relative_motion(:, 1))
grid on
axis equal
xlabel('N-axis')
ylabel('R-axis')

subplot(2,2,3)
plot(relative_motion(:, 2), relative_motion(:, 3))
grid on
axis equal
xlabel('T-axis')
ylabel('N-axis')

subplot(2,2,4)
plot3(relative_motion(:, 1), relative_motion(:, 2), relative_motion(:, 3))
grid on
axis equal
view(3)
xlabel('R-axis')
ylabel('T-axis')
zlabel('N-axis')

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
        e_d * sin(omega_d) - e_c * cos(omega_c);
        inc_d - inc_c;
        (RAAN_d - RAAN_c) * sin(inc_c)]

%% Part f: Propagating ROEs using linear mapping



%% Functions

function M = STM_YA(f, e, t, n, M_1)
    k = 1 + e * cos(f);
    k_prime = - e * sin(f);
    eta = sqrt(1 - e^2);
    tau = n * t / eta^3;
    
    M_2 = [      1/k + 3/2*k_prime*tau,          sin(f),            cos(f),       0,          0,          0;
                            -3/2*k*tau,  (1+1/k)*cos(f),   -(1+1/k)*sin(f),     1/k,          0,          0;
                                     0,               0,                 0,       0, 1/k*sin(f), 1/k*cos(f);
           k_prime/2-3/2*k^2*(k-1)*tau,      k^2*cos(f),       -k^2*sin(f),       0,          0,          0;
               -3/2*(k+k^2*k_prime*tau), -(k^2+1)*sin(f), -e-(k^2+1)*cos(f), -k_prime,          0,          0;
                                     0,               0,                 0,       0,   e+cos(f),    -sin(f)];
    M = M_1 * M_2;
end