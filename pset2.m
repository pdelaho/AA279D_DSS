%% Problem 1  

close all; clc; clear;
path_config;
mu = 398600.435436; % km^3/s^2
R_E = 6378; % km

%% Part a: Chief and deputy orbital elements

% For the chief
% Since there are no perturbations, mean elements = osculating elements
a_c = 36943; % km, semi-major axis
% e_c = 0.8111; % eccentricity
e_c = 0.0001;
inc_c = deg2rad(59); % inclination
omega_c = deg2rad(188); % argument of periapsis
RAAN_c = deg2rad(84); % right ascension of the ascending node
nu_c = deg2rad(0); % true anomaly
T = 2 * pi * sqrt(a_c^3 / mu);
n_c = sqrt(mu / a_c^3);

% % For the deputy
% delta_a = 0;
% delta_e = - e_c * 0.01 * a_c;
% delta_i = deg2rad(0.5) * a_c;
% delta_omega = 0;
% delta_RAAN = - deg2rad(1) * a_c;
% delta_nu = deg2rad(-0.016) * a_c;

% For the deputy
delta_a = 0;
delta_e = - e_c * 0.1 * a_c;
delta_i = deg2rad(0.0005) * a_c;
delta_omega = 0;
delta_RAAN = - deg2rad(0.0001) * a_c;
delta_nu = deg2rad(-0.0016) * a_c;

a_d = a_c + delta_a;
e_d = e_c + delta_e / a_c;
inc_d = inc_c + delta_i / a_c;
omega_d = omega_c + delta_omega / a_c;
RAAN_d = RAAN_c + delta_RAAN / a_c;
nu_d = nu_c + delta_nu / a_c;

%% Part b: Numerical integration of NERM

[pos_c, vel_c] = OEtoECI(a_c, e_c, inc_c, omega_c, RAAN_c, nu_c, mu);
[pos_d, vel_d] = OEtoECI(a_d, e_d, inc_d, omega_d, RAAN_d, nu_d, mu);

R_vec = pos_c / norm(pos_c);
N_vec = cross(pos_c, vel_c) / norm(cross(pos_c, vel_c));
T_vec = cross(N_vec, R_vec);
rot = [R_vec'; T_vec'; N_vec'];
theta_dot = sqrt(mu / (a_c^3 * (1 - e_c^2)^3)) * (1 + e_c * cos(nu_c))^2;
theta_dot_vec = [0 0 theta_dot];
r0_dot = norm(pos_c) * e_c * sin(nu_c) / (1 + e_c * cos(nu_c)) * sqrt(mu / (a_c^3 * (1 - e_c^2)^3)) * (1 + e_c * cos(nu_c))^2;

rel_pos = pos_d - pos_c;
rel_pos_RTN = (rot * rel_pos)';
rel_vel = vel_d - vel_c;
rel_vel_RTN = (rot * rel_vel)' - cross(theta_dot_vec, rel_pos_RTN);
% vel_RTN_c = (rot * vel_c)' - cross(theta_dot_vec, rot * pos_c);

% init_cond_rel = [rel_pos_RTN rel_vel_RTN norm(pos_c) norm(vel_RTN_c) nu_c theta_dot];
init_cond_rel = [rel_pos_RTN rel_vel_RTN norm(pos_c) r0_dot nu_c theta_dot];
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-15, 'MaxStep', 100);
N = 10000;
tspan = linspace(0, 10 * T, N);
[t, y] = ode89(@(t, state) NERM(t, state, mu), tspan, init_cond_rel, options);

% Plots for the position
figure
subplot(6,4,[1,2,5,6])
plot(t / T, y(:, 1) / a_c)
grid on
sgtitle("Normalized relative position of the deputy in the RTN frame based on the chief using NERM (vs time and 2D/3D plots)")
ylabel('x/a_0 [-]')

subplot(6,4,[9,10,13,14])
plot(t / T, y(:, 2) / a_c)
grid on
ylabel('y/a_0 [-]')

subplot(6,4,[17,18,21,22])
plot(t / T, y(:, 3) / a_c)
grid on
ylabel('z/a_0 [-]')
xlabel("Number of chief's orbit")

subplot(6,4,[3,7,11])
plot(y(:, 1) / a_c, y(:, 2) / a_c)
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')

subplot(6,4,[4,8,12])
plot(y(:, 1) / a_c, y(:, 3) / a_c)
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[15,19,23])
plot(y(:, 2) / a_c, y(:, 3) / a_c)
grid on
axis square
xlabel('y/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[16,20,24])
plot3(y(:, 1) / a_c, y(:, 2) / a_c, y(:, 3) / a_c)
axis square
grid on
view(3)
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')
zlabel('z/a_0 [-]')

% Plots for the velocity
figure
subplot(6,4,[1,2,5,6])
plot(t / T, y(:, 4) / (a_c * n_c))
grid on
sgtitle("Normalized relative velocity of the deputy in the RTN frame based on the chief using NERM (vs time and 2D/3D plots)")
ylabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[9,10,13,14])
plot(t / T, y(:, 5) / (a_c * n_c))
grid on
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[17,18,21,22])
plot(t / T, y(:, 6) / (a_c * n_c))
grid on
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')
xlabel("Number of chief's orbit")


subplot(6,4,[3,7,11])
plot(y(:, 4) / (a_c * n_c), y(:, 5) / (a_c * n_c))
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[4,8,12])
plot(y(:, 4) / (a_c * n_c), y(:, 6) / (a_c * n_c))
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[15,19,23])
plot(y(:, 5) / (a_c * n_c), y(:, 6) / (a_c * n_c))
grid on
axis square
xlabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[16,20,24])
plot3(y(:, 4) / (a_c * n_c), y(:, 5) / (a_c * n_c), y(:, 6) / (a_c * n_c))
axis square
grid on
view(3)
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
zlabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

%% Part c: Numerical integration of FODE

init_cond_abs = [pos_c', vel_c', pos_d', vel_d'];
[t_abs, y_abs] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, init_cond_abs, options);

fode_RTN = zeros(N, 6);

for i=1:N
    pos_c_temp = y_abs(i, 1:3);
    vel_c_temp = y_abs(i, 4:6);
    rho_temp = y_abs(i, 7:9) - pos_c_temp;
    rho_dot_temp = y_abs(i, 10:12) - vel_c_temp;
    [a, e, inc, omega, RAAN, nu] = Keplerian_elements(y_abs(i, 1:6), mu);

    % Computing the RTN frame centered on the chief spacecraft
    R_vec_temp = pos_c_temp / norm(pos_c_temp);
    N_vec_temp = cross(pos_c_temp, vel_c_temp) / norm(cross(pos_c_temp, vel_c_temp));
    T_vec_temp = cross(N_vec_temp, R_vec_temp);
    rot_temp = [R_vec_temp; T_vec_temp; N_vec_temp];
    theta_dot_temp = [0, 0, - sqrt(mu / (a^3 * (1 - e^2)^3)) * (1 + e * cos(nu))^2];
    
    % Transforming the position and velocities in the RTN frame
    rho_RTN = rot_temp * rho_temp';
    rho_dot_RTN = rot_temp * rho_dot_temp' + cross(theta_dot_temp, rho_RTN)';

    fode_RTN(i, 1:3) = rho_RTN;
    fode_RTN(i, 4:6) = rho_dot_RTN;
end

% Plots for the position
figure
subplot(6,4,[1,2,5,6])
hold on
plot(t / T, y(:, 1) / a_c, 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN(:, 1) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
legend('NERM', 'FODE')
sgtitle("Normalized relative position of the deputy in the RTN frame based on the chief using FODE (vs time and 2D/3D plots)")
ylabel('x/a_0 [-]')

subplot(6,4,[9,10,13,14])
hold on
plot(t / T, y(:, 2) / a_c, 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN(:, 2) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
ylabel('y/a_0 [-]')

subplot(6,4,[17,18,21,22])
hold on
plot(t / T, y(:, 3) / a_c, 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN(:, 3) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
ylabel('z/a_0 [-]')
xlabel("Number of chief's orbit")

subplot(6,4,[3,7,11])
hold on
plot(y(:, 1) / a_c, y(:, 2) / a_c, 'b-', 'LineWidth', 1.5)
plot(fode_RTN(:, 1) / a_c, fode_RTN(:, 2) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')

subplot(6,4,[4,8,12])
hold on
plot(y(:, 1) / a_c, y(:, 3) / a_c, 'b-', 'LineWidth', 1.5)
plot(fode_RTN(:, 1) / a_c, fode_RTN(:, 3) / a_c, 'r--', 'LineWidth', 1.5)
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[15,19,23])
hold on
plot(y(:, 2) / a_c, y(:, 3) / a_c, 'b-', 'LineWidth', 1.5)
plot(fode_RTN(:, 2) / a_c, fode_RTN(:, 3) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('y/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[16,20,24])
hold on
plot3(y(:, 1) / a_c, y(:, 2) / a_c, y(:, 3) / a_c, 'b-', 'LineWidth', 1.5)
plot3(fode_RTN(:, 1) / a_c, fode_RTN(:, 2) / a_c, fode_RTN(:, 3) / a_c, 'r--', 'LineWidth', 1.5)
hold off
axis square
grid on
view(3)
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')
zlabel('z/a_0 [-]')

% Plots for the velocity
figure
subplot(6,4,[1,2,5,6])
hold on
plot(t / T, y(:, 4) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN(:, 4) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
legend('NERM', 'FODE')
sgtitle("Normalized relative velocity of the deputy in the RTN frame based on the chief using FODE (vs time and 2D/3D plots)")
ylabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[9,10,13,14])
hold on
plot(t / T, y(:, 5) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN(:, 5) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[17,18,21,22])
hold on
plot(t / T, y(:, 6) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN(:, 6) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')
xlabel("Number of chief's orbit")

subplot(6,4,[3,7,11])
hold on
plot(y(:, 4) / (a_c * n_c), y(:, 5) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(fode_RTN(:, 4) / (a_c * n_c), fode_RTN(:, 5) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[4,8,12])
hold on
plot(y(:, 4) / (a_c * n_c), y(:, 6) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(fode_RTN(:, 4) / (a_c * n_c), fode_RTN(:, 6) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[15,19,23])
hold on
plot(y(:, 5) / (a_c * n_c), y(:, 6) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(fode_RTN(:, 5) / (a_c * n_c), fode_RTN(:, 6) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[16,20,24])
hold on
plot3(y(:, 4) / (a_c * n_c), y(:, 5) / (a_c * n_c), y(:, 6) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot3(fode_RTN(:, 4) / (a_c * n_c), fode_RTN(:, 5) / (a_c * n_c), fode_RTN(:, 6) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
axis square
grid on
view(3)
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
zlabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

% Modeling the Earth
[X, Y, Z] = sphere(10);

figure
surf(X * R_E, Y * R_E, Z * R_E)
hold on
plot3(y_abs(:, 1), y_abs(:, 2), y_abs(:, 3))
plot3(y_abs(:, 7), y_abs(:, 8), y_abs(:, 9))
axis equal
grid on
hold off
view(3)
legend('Earth', 'Chief orbit', 'Deputy orbit')
xlabel('X-axis [km]')
ylabel('Y-axis [km]')
zlabel('Z-axis [km]')
title("Chief's and deputy's orbits around the Earth in ECI frame")

%% Part d: Numerical errors + new initial conditions

% Checking it's only numerical errors between parts b and c
figure
sgtitle('Error in relative position and velocity in RTN frame between NERM and FODE')
subplot(3,2,1)
plot(t_abs / T, abs(y(:, 1) - fode_RTN(:, 1)))
grid on
ylabel('R-component of position [km]')

subplot(3,2,3)
plot(t_abs / T, abs(y(:, 2) - fode_RTN(:, 2)))
grid on
ylabel('T-component of position [km]')

subplot(3,2,5)
plot(t_abs / T, abs(y(:, 3) - fode_RTN(:, 3)))
grid on
ylabel('N-component of position [km]')
xlabel("Number of chief's orbit")

subplot(3,2,2)
plot(t_abs / T, abs(y(:, 4) - fode_RTN(:, 4)))
grid on
ylabel('R-component of velocity [km/s]')

subplot(3,2,4)
plot(t_abs / T, abs(y(:, 5) - fode_RTN(:, 5)))
grid on
ylabel('T-component of velocity [km/s]')

subplot(3,2,6)
plot(t_abs / T, abs(y(:, 6) - fode_RTN(:, 6)))
grid on
ylabel('N-component of velocity [km/s]')
xlabel("Number of chief's orbit")

% New set of initial conditions with a difference in semi-major axis
delta_a2 = 10;
delta_e2 = 0;
delta_i2 = 0;
delta_omega2 = deg2rad(7) * a_c;
delta_RAAN2 = deg2rad(5) * a_c;
delta_nu2 = deg2rad(-0.016) * a_c;
% delta_nu2 = 0;

a_d2 = a_c + delta_a2;
e_d2 = e_c + delta_e2 / a_c;
inc_d2 = inc_c + delta_i2 / a_c;
omega_d2 = omega_c + delta_omega2 / a_c;
RAAN_d2 = RAAN_c + delta_RAAN2 / a_c;
nu_d2 = nu_c + delta_nu2 / a_c;

[pos_d2, vel_d2] = OEtoECI(a_d2, e_d2, inc_d2, omega_d2, RAAN_d2, nu_d2, mu);

rel_pos2 = pos_d2 - pos_c;
rel_pos_RTN2 = (rot * rel_pos2)';
rel_vel2 = vel_d2 - vel_c;
rel_vel_RTN2 = (rot * rel_vel2)' - cross(theta_dot_vec, rel_pos_RTN2);
% vel_RTN_c2 = (rot * vel_c)' - cross(theta_dot_vec, rot * pos_c);

% init_cond_rel = [rel_pos_RTN rel_vel_RTN norm(pos_c) norm(vel_RTN_c) nu_c theta_dot];
init_cond_rel2 = [rel_pos_RTN2 rel_vel_RTN2 norm(pos_c) r0_dot nu_c theta_dot];
[t2, y2] = ode89(@(t, state) NERM(t, state, mu), tspan, init_cond_rel2, options);

init_cond_abs2 = [pos_c', vel_c', pos_d2', vel_d2'];
[t_abs2, y_abs2] = ode89(@(t, state) FODE_2sats(t, state, mu), tspan, init_cond_abs2, options);

fode_RTN2 = zeros(N, 6);

for i=1:N
    pos_c_temp = y_abs2(i, 1:3);
    vel_c_temp = y_abs2(i, 4:6);
    rho_temp = y_abs2(i, 7:9) - pos_c_temp;
    rho_dot_temp = y_abs2(i, 10:12) - vel_c_temp;
    [a, e, inc, omega, RAAN, nu] = Keplerian_elements(y_abs2(i, 1:6), mu);

    % Computing the RTN frame centered on the chief spacecraft
    R_vec_temp = pos_c_temp / norm(pos_c_temp);
    N_vec_temp = cross(pos_c_temp, vel_c_temp) / norm(cross(pos_c_temp, vel_c_temp));
    T_vec_temp = cross(N_vec_temp, R_vec_temp);
    rot_temp = [R_vec_temp; T_vec_temp; N_vec_temp];
    theta_dot_temp = [0, 0, - sqrt(mu / (a^3 * (1 - e^2)^3)) * (1 + e * cos(nu))^2];
    
    % Transforming the position and velocities in the RTN frame
    rho_RTN = rot_temp * rho_temp';
    rho_dot_RTN = rot_temp * rho_dot_temp' + cross(theta_dot_temp, rho_RTN)';

    fode_RTN2(i, 1:3) = rho_RTN;
    fode_RTN2(i, 4:6) = rho_dot_RTN;
end

% Plots for the position
figure
subplot(6,4,[1,2,5,6])
hold on
plot(t / T, y2(:, 1) / a_c, 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN2(:, 1) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
legend('NERM', 'FODE')
sgtitle("Normalized relative position of the deputy in the RTN frame based on the chief using FODE (vs time and 2D/3D plots)")
ylabel('x/a_0 [-]')

subplot(6,4,[9,10,13,14])
hold on
plot(t / T, y2(:, 2) / a_c, 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN2(:, 2) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
ylabel('y/a_0 [-]')

subplot(6,4,[17,18,21,22])
hold on
plot(t / T, y2(:, 3) / a_c, 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN2(:, 3) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
ylabel('z/a_0 [-]')
xlabel("Chief's orbital periods")

subplot(6,4,[3,7,11])
hold on
plot(y2(:, 1) / a_c, y2(:, 2) / a_c, 'b-', 'LineWidth', 1.5)
plot(fode_RTN2(:, 1) / a_c, fode_RTN2(:, 2) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')

subplot(6,4,[4,8,12])
hold on
plot(y2(:, 1) / a_c, y2(:, 3) / a_c, 'b-', 'LineWidth', 1.5)
plot(fode_RTN2(:, 1) / a_c, fode_RTN2(:, 3) / a_c, 'r--', 'LineWidth', 1.5)
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[15,19,23])
hold on
plot(y2(:, 2) / a_c, y2(:, 3) / a_c, 'b-', 'LineWidth', 1.5)
plot(fode_RTN2(:, 2) / a_c, fode_RTN2(:, 3) / a_c, 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('y/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[16,20,24])
hold on
plot3(y2(:, 1) / a_c, y2(:, 2) / a_c, y2(:, 3) / a_c, 'b-', 'LineWidth', 1.5)
plot3(fode_RTN2(:, 1) / a_c, fode_RTN2(:, 2) / a_c, fode_RTN2(:, 3) / a_c, 'r--', 'LineWidth', 1.5)
hold off
axis square
grid on
view(3)
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')
zlabel('z/a_0 [-]')

% Plots for the velocity
figure
subplot(6,4,[1,2,5,6])
hold on
plot(t / T, y2(:, 4) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN2(:, 4) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
legend('NERM', 'FODE')
sgtitle("Normalized relative velocity of the deputy in the RTN frame based on the chief using FODE (vs time and 2D/3D plots)")
ylabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[9,10,13,14])
hold on
plot(t / T, y2(:, 5) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN2(:, 5) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[17,18,21,22])
hold on
plot(t / T, y2(:, 6) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(t / T, fode_RTN2(:, 6) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')
xlabel("Chief's orbital periods")

subplot(6,4,[3,7,11])
hold on
plot(y2(:, 4) / (a_c * n_c), y2(:, 5) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(fode_RTN2(:, 4) / (a_c * n_c), fode_RTN2(:, 5) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[4,8,12])
hold on
plot(y2(:, 4) / (a_c * n_c), y2(:, 6) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(fode_RTN2(:, 4) / (a_c * n_c), fode_RTN2(:, 6) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[15,19,23])
hold on
plot(y2(:, 5) / (a_c * n_c), y2(:, 6) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot(fode_RTN2(:, 5) / (a_c * n_c), fode_RTN2(:, 6) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
grid on
axis square
xlabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[16,20,24])
hold on
plot3(y2(:, 4) / (a_c * n_c), y2(:, 5) / (a_c * n_c), y2(:, 6) / (a_c * n_c), 'b-', 'LineWidth', 1.5)
plot3(fode_RTN2(:, 4) / (a_c * n_c), fode_RTN2(:, 5) / (a_c * n_c), fode_RTN2(:, 6) / (a_c * n_c), 'r--', 'LineWidth', 1.5)
hold off
axis square
grid on
view(3)
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
zlabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

% Modeling the Earth
[X, Y, Z] = sphere(10);

figure
surf(X * R_E, Y * R_E, Z * R_E)
hold on
plot3(y_abs2(:, 1), y_abs2(:, 2), y_abs2(:, 3))
plot3(y_abs2(:, 7), y_abs2(:, 8), y_abs2(:, 9))
axis equal
grid on
hold off
view(3)
legend('Earth', 'Chief orbit', 'Deputy orbit')
xlabel('X-axis [km]')
ylabel('Y-axis [km]')
zlabel('Z-axis [km]')
title("Chief's and deputy's orbits around the Earth in ECI frame")

% Checking it's only numerical errors between parts b and c
figure
sgtitle('Error in relative position and velocity in RTN frame between NERM and FODE')
subplot(3,2,1)
plot(t_abs2 / T, abs(y2(:, 1) - fode_RTN2(:, 1)))
grid on
ylabel('R-component of position [km]')

subplot(3,2,3)
plot(t_abs2 / T, abs(y2(:, 2) - fode_RTN2(:, 2)))
grid on
ylabel('T-component of position [km]')

subplot(3,2,5)
plot(t_abs2 / T, abs(y2(:, 3) - fode_RTN2(:, 3)))
grid on
ylabel('N-component of position [km]')
xlabel("Number of chief's orbit")

subplot(3,2,2)
plot(t_abs2 / T, abs(y2(:, 4) - fode_RTN2(:, 4)))
grid on
ylabel('R-component of velocity [km/s]')

subplot(3,2,4)
plot(t_abs2 / T, abs(y2(:, 5) - fode_RTN2(:, 5)))
grid on
ylabel('T-component of velocity [km/s]')

subplot(3,2,6)
plot(t_abs2 / T, abs(y2(:, 6) - fode_RTN2(:, 6)))
grid on
ylabel('N-component of velocity [km/s]')
xlabel("Number of chief's orbit")

%% Part f: Maneuver to establish a bounded periodic relative motion

E_d2 = 2 * atan(sqrt((1 - e_d2) / (1 + e_d2)) * tan(nu_d2 / 2));
M_d2 = E_d2 - e_d2 * sin(E_d2);
n_d2 = sqrt(mu / a_d2^3);
T_d2 = 2 * pi / n_d2;
delta_t = (0 - M_d2) / n_d2; % flight time between initial mean anomaly and deputy's periapsis
tspan_before = linspace(0, 3 * T_d2 + delta_t, 3000);
[t_man_bef, y_man_bef] = ode89(@(t, state) NERM(t, state, mu), tspan_before, init_cond_rel2, options);

delta_v = - delta_a2 * n_d2 * sqrt(1 - e_d2) / (2 * sqrt(1 + e_d2));
init_cond_maneuver = y_man_bef(end,:) + [0 0 0 0 delta_v 0 0 0 0 0];

tspan_after = linspace(tspan_before(end), 10 * T, 7000);
[t_man_aft, y_man_aft] = ode89(@(t, state) NERM(t, state, mu), tspan_after, init_cond_maneuver, options);

init_cond_bef_maneuver = [init_cond_rel2 init_cond_abs2];
[t_man_bef_2, y_man_bef_2] = ode89(@(t, state) NERM_maneuver(t, state, mu), tspan_before, init_cond_bef_maneuver, options);

maneuver_RTN_d = [0 delta_v 0];
[a, e, inc, omega, RAAN, nu] = Keplerian_elements(y_man_bef_2(end, 11:16), mu);
% Computing the RTN frame centered on the chief spacecraft
R_vec_c = y_man_bef_2(end, 11:13) / norm(y_man_bef_2(end, 11:13));
N_vec_c = cross(y_man_bef_2(end, 11:13), y_man_bef_2(end, 14:16)) / norm(cross(y_man_bef_2(end, 11:13), y_man_bef_2(end, 14:16)));
T_vec_c = cross(N_vec_c, R_vec_c);
rot_ECItoRTN_c = [R_vec_c; T_vec_c; N_vec_c];
theta_dot_c = [0, 0, sqrt(mu / (a^3 * (1 - e^2)^3)) * (1 + e * cos(nu))^2];

[a, e, inc, omega, RAAN, nu] = Keplerian_elements(y_man_bef_2(end, 17:22), mu);
% Computing the RTN frame centered on the chief spacecraft
R_vec_d = y_man_bef_2(end, 17:19) / norm(y_man_bef_2(end, 17:19));
N_vec_d = cross(y_man_bef_2(end, 17:19), y_man_bef_2(end, 20:22)) / norm(cross(y_man_bef_2(end, 17:19), y_man_bef_2(end, 20:22)));
T_vec_d = cross(N_vec_d, R_vec_d);
rot_ECItoRTN_d = [R_vec_d; T_vec_d; N_vec_d];
% theta_dot_d = [0, 0, - sqrt(mu / (a^3 * (1 - e^2)^3)) * (1 + e * cos(nu))^2];

maneuver_RTN_c = rot_ECItoRTN_c * rot_ECItoRTN_d' * maneuver_RTN_d';

init_cond_aft_maneuver = y_man_bef_2(end, 1:10) + [0 0 0 maneuver_RTN_c' 0 0 0 0];
[t_man_aft_2, y_man_aft_2] = ode89(@(t, state) NERM(t, state, mu), tspan_after, init_cond_aft_maneuver, options);

[a, e, inc, omega, RAAN, nu] = Keplerian_elements(y_abs2(end, 1:6), mu);
% Computing the RTN frame centered on the chief spacecraft
R_vec_c = y_abs2(end, 1:3) / norm(y_abs2(end, 1:3));
N_vec_c = cross(y_abs2(end, 1:3), y_abs2(end, 4:6)) / norm(cross(y_abs2(end, 1:3), y_abs2(end, 4:6)));
T_vec_c = cross(N_vec_c, R_vec_c);
rot_ECItoRTN_c = [R_vec_c; T_vec_c; N_vec_c];
theta_dot_c = [0, 0, sqrt(mu / (a^3 * (1 - e^2)^3)) * (1 + e * cos(nu))^2];

rel_pos_ECI = rot_ECItoRTN_c' * y_man_aft(end, 1:3)';
rel_vel_ECI = rot_ECItoRTN_c' * (y_man_aft(end, 4:6)' + cross(theta_dot_c, y_man_aft(end, 1:3))');
pos_d_ECI = y_abs2(end, 1:3) + rel_pos_ECI';
vel_d_ECI = y_abs2(end, 4:6) + rel_vel_ECI';
[a, e, inc, omega, RAAN, nu] = Keplerian_elements([pos_d_ECI vel_d_ECI], mu);

a - a_c
e - e_c
delta_e_man = 2 * sqrt(1 - e_d2^2) / (n_d2 * a_d2) * delta_v

rel_pos_ECI = rot_ECItoRTN_c' * y_man_aft_2(end, 1:3)';
rel_vel_ECI = rot_ECItoRTN_c' * (y_man_aft_2(end, 4:6)' + cross(theta_dot_c, y_man_aft_2(end, 1:3))');
pos_d_ECI = y_abs2(end, 1:3) + rel_pos_ECI';
vel_d_ECI = y_abs2(end, 4:6) + rel_vel_ECI';
[a, e, inc, omega, RAAN, nu] = Keplerian_elements([pos_d_ECI vel_d_ECI], mu);

a - a_c
e - e_c

% Plots for the position
figure
subplot(6,4,[1,2,5,6])
hold on
plot(t_man_bef / T, y_man_bef(:, 1) / a_c)
plot(t_man_aft / T, y_man_aft(:, 1) / a_c)
hold off
legend('Before \DeltaV', 'After \DeltaV')
grid on
sgtitle("Normalized relative position of the deputy in the RTN frame based on the chief using NERM (vs time and 2D/3D plots)")
ylabel('x/a_0 [-]')

subplot(6,4,[9,10,13,14])
hold on
plot(t_man_bef / T, y_man_bef(:, 2) / a_c)
plot(t_man_aft / T, y_man_aft(:, 2) / a_c)
hold off
grid on
ylabel('y/a_0 [-]')

subplot(6,4,[17,18,21,22])
hold on
plot(t_man_bef / T, y_man_bef(:, 3) / a_c)
plot(t_man_aft / T, y_man_aft(:, 3) / a_c)
hold off
grid on
ylabel('z/a_0 [-]')
xlabel("Number of chief's orbit")

subplot(6,4,[3,7,11])
hold on
plot(y_man_bef(:, 1) / a_c, y_man_bef(:, 2) / a_c)
plot(y_man_aft(:, 1) / a_c, y_man_aft(:, 2) / a_c)
hold off
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')

subplot(6,4,[4,8,12])
hold on
plot(y_man_bef(:, 1) / a_c, y_man_bef(:, 3) / a_c)
plot(y_man_aft(:, 1) / a_c, y_man_aft(:, 3) / a_c)
hold off
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[15,19,23])
hold on
plot(y_man_bef(:, 2) / a_c, y_man_bef(:, 3) / a_c)
plot(y_man_aft(:, 2) / a_c, y_man_aft(:, 3) / a_c)
hold off
grid on
axis square
xlabel('y/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[16,20,24])
hold on
plot3(y_man_bef(:, 1) / a_c, y_man_bef(:, 2) / a_c, y_man_bef(:, 3) / a_c)
plot3(y_man_aft(:, 1) / a_c, y_man_aft(:, 2) / a_c, y_man_aft(:, 3) / a_c)
hold off
axis equal
grid on
view(3)
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')
zlabel('z/a_0 [-]')

% Plots for the velocity
figure
subplot(6,4,[1,2,5,6])
hold on
plot(t_man_bef / T, y_man_bef(:, 4) / (a_c * n_c))
plot(t_man_aft / T, y_man_aft(:, 4) / (a_c * n_c))
hold off
legend('Before \DeltaV', 'After \DeltaV')
grid on
sgtitle("Normalized relative velocity of the deputy in the RTN frame based on the chief using NERM (vs time and 2D/3D plots)")
ylabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[9,10,13,14])
hold on
plot(t_man_bef / T, y_man_bef(:, 5) / (a_c * n_c))
plot(t_man_aft / T, y_man_aft(:, 5) / (a_c * n_c))
hold off
grid on
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[17,18,21,22])
hold on
plot(t_man_bef / T, y_man_bef(:, 6) / (a_c * n_c))
plot(t_man_aft / T, y_man_aft(:, 6) / (a_c * n_c))
hold off
grid on
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')
xlabel("Number of chief's orbit")


subplot(6,4,[3,7,11])
hold on
plot(y_man_bef(:, 4) / (a_c * n_c), y_man_bef(:, 5) / (a_c * n_c))
plot(y_man_aft(:, 4) / (a_c * n_c), y_man_aft(:, 5) / (a_c * n_c))
hold off
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[4,8,12])
hold on
plot(y_man_bef(:, 4) / (a_c * n_c), y_man_bef(:, 6) / (a_c * n_c))
plot(y_man_aft(:, 4) / (a_c * n_c), y_man_aft(:, 6) / (a_c * n_c))
hold off
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[15,19,23])
hold on
plot(y_man_bef(:, 5) / (a_c * n_c), y_man_bef(:, 6) / (a_c * n_c))
plot(y_man_aft(:, 5) / (a_c * n_c), y_man_aft(:, 6) / (a_c * n_c))
hold off
grid on
axis square
xlabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[16,20,24])
hold on
plot3(y_man_bef(:, 4) / (a_c * n_c), y_man_bef(:, 5) / (a_c * n_c), y_man_bef(:, 6) / (a_c * n_c))
plot3(y_man_aft(:, 4) / (a_c * n_c), y_man_aft(:, 5) / (a_c * n_c), y_man_aft(:, 6) / (a_c * n_c))
hold off
axis equal
grid on
view(3)
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
zlabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

% Zoom
figure
hold on
plot(y_man_bef(:, 4) / (a_c * n_c), y_man_bef(:, 5) / (a_c * n_c))
plot(y_man_aft(:, 4) / (a_c * n_c), y_man_aft(:, 5) / (a_c * n_c))
hold off
grid on
axis square
xlim([-0.18, -0.15])
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

% Plots for the position
figure
subplot(6,4,[1,2,5,6])
hold on
plot(t_man_bef_2 / T, y_man_bef_2(:, 1) / a_c)
plot(t_man_aft_2 / T, y_man_aft_2(:, 1) / a_c)
hold off
legend('Before \DeltaV', 'After \DeltaV')
grid on
sgtitle("Normalized relative position of the deputy in the RTN frame based on the chief using NERM (vs time and 2D/3D plots)")
ylabel('x/a_0 [-]')

subplot(6,4,[9,10,13,14])
hold on
plot(t_man_bef_2 / T, y_man_bef_2(:, 2) / a_c)
plot(t_man_aft_2 / T, y_man_aft_2(:, 2) / a_c)
hold off
grid on
ylabel('y/a_0 [-]')

subplot(6,4,[17,18,21,22])
hold on
plot(t_man_bef_2 / T, y_man_bef_2(:, 3) / a_c)
plot(t_man_aft_2 / T, y_man_aft_2(:, 3) / a_c)
hold off
grid on
ylabel('z/a_0 [-]')
xlabel("Number of chief's orbit")

subplot(6,4,[3,7,11])
hold on
plot(y_man_bef_2(:, 1) / a_c, y_man_bef_2(:, 2) / a_c)
plot(y_man_aft_2(:, 1) / a_c, y_man_aft_2(:, 2) / a_c)
hold off
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')

subplot(6,4,[4,8,12])
hold on
plot(y_man_bef_2(:, 1) / a_c, y_man_bef_2(:, 3) / a_c)
plot(y_man_aft_2(:, 1) / a_c, y_man_aft_2(:, 3) / a_c)
hold off
grid on
axis square
xlabel('x/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[15,19,23])
hold on
plot(y_man_bef_2(:, 2) / a_c, y_man_bef_2(:, 3) / a_c)
plot(y_man_aft_2(:, 2) / a_c, y_man_aft_2(:, 3) / a_c)
hold off
grid on
axis square
xlabel('y/a_0 [-]')
ylabel('z/a_0 [-]')

subplot(6,4,[16,20,24])
hold on
plot3(y_man_bef_2(:, 1) / a_c, y_man_bef_2(:, 2) / a_c, y_man_bef_2(:, 3) / a_c)
plot3(y_man_aft_2(:, 1) / a_c, y_man_aft_2(:, 2) / a_c, y_man_aft_2(:, 3) / a_c)
hold off
axis equal
grid on
view(3)
xlabel('x/a_0 [-]')
ylabel('y/a_0 [-]')
zlabel('z/a_0 [-]')

% Plots for the velocity
figure
subplot(6,4,[1,2,5,6])
hold on
plot(t_man_bef_2 / T, y_man_bef_2(:, 4) / (a_c * n_c))
plot(t_man_aft_2 / T, y_man_aft_2(:, 4) / (a_c * n_c))
hold off
legend('Before \DeltaV', 'After \DeltaV')
grid on
sgtitle("Normalized relative velocity of the deputy in the RTN frame based on the chief using NERM (vs time and 2D/3D plots)")
ylabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[9,10,13,14])
hold on
plot(t_man_bef_2 / T, y_man_bef_2(:, 5) / (a_c * n_c))
plot(t_man_aft_2 / T, y_man_aft_2(:, 5) / (a_c * n_c))
hold off
grid on
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[17,18,21,22])
hold on
plot(t_man_bef_2 / T, y_man_bef_2(:, 6) / (a_c * n_c))
plot(t_man_aft_2 / T, y_man_aft_2(:, 6) / (a_c * n_c))
hold off
grid on
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')
xlabel("Number of chief's orbit")


subplot(6,4,[3,7,11])
hold on
plot(y_man_bef_2(:, 4) / (a_c * n_c), y_man_bef_2(:, 5) / (a_c * n_c))
plot(y_man_aft_2(:, 4) / (a_c * n_c), y_man_aft_2(:, 5) / (a_c * n_c))
hold off
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[4,8,12])
hold on
plot(y_man_bef_2(:, 4) / (a_c * n_c), y_man_bef_2(:, 6) / (a_c * n_c))
plot(y_man_aft_2(:, 4) / (a_c * n_c), y_man_aft_2(:, 6) / (a_c * n_c))
hold off
grid on
axis square
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[15,19,23])
hold on
plot(y_man_bef_2(:, 5) / (a_c * n_c), y_man_bef_2(:, 6) / (a_c * n_c))
plot(y_man_aft_2(:, 5) / (a_c * n_c), y_man_aft_2(:, 6) / (a_c * n_c))
hold off
grid on
axis square
xlabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

subplot(6,4,[16,20,24])
hold on
plot3(y_man_bef_2(:, 4) / (a_c * n_c), y_man_bef_2(:, 5) / (a_c * n_c), y_man_bef_2(:, 6) / (a_c * n_c))
plot3(y_man_aft_2(:, 4) / (a_c * n_c), y_man_aft_2(:, 5) / (a_c * n_c), y_man_aft_2(:, 6) / (a_c * n_c))
hold off
axis equal
grid on
view(3)
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')
zlabel('$\dot{z}/(a_0n_0) [-]$', 'Interpreter', 'latex')

% Zoom
figure
hold on
plot(y_man_bef_2(:, 4) / (a_c * n_c), y_man_bef_2(:, 5) / (a_c * n_c))
plot(y_man_aft_2(:, 4) / (a_c * n_c), y_man_aft_2(:, 5) / (a_c * n_c))
hold off
grid on
axis square
xlim([-0.18, -0.15])
xlabel('$\dot{x}/(a_0n_0) [-]$', 'Interpreter', 'latex')
ylabel('$\dot{y}/(a_0n_0) [-]$', 'Interpreter', 'latex')

%% Functions

function [pos_inertial, vel_inertial] = OEtoECI(a, e, inc, omega, RAAN, true_anom, mu)
    cosE = (e + cos(true_anom)) / (1 + e * cos(true_anom));
    sinE = sin(true_anom) * sqrt(1 - e^2) / (1 + e * cos(true_anom));
    n = sqrt(mu / a^3);
    
    % Computing the position and velocity vectors in the perifocal frame
    pos_perifocal = [a * (cosE - e) a * sqrt(1 - e^2) * sinE 0];
    vel_perifocal = a * n / (1 - e * cosE) * [-sinE sqrt(1 - e^2)*cosE 0];
    
    % Computing the rotation matrix between perifocal and inertial frame
    rotRAAN = [cos(RAAN) sin(-RAAN) 0;
               -sin(-RAAN) cos(RAAN) 0;
               0 0 1];
    roti = [1 0 0;
            0 cos(inc) sin(-inc);
            0 -sin(-inc) cos(inc)];
    rotomega = [cos(omega) sin(-omega) 0;
                -sin(-omega) cos(omega) 0;
                0 0 1];
    rot_perifocalTOinertial = rotRAAN * roti * rotomega;

    % Rotating the position and velocity vectors
    pos_inertial = rot_perifocalTOinertial * pos_perifocal';
    vel_inertial = rot_perifocalTOinertial * vel_perifocal';
end

function statedot = NERM(t, state, mu)
    statedot = zeros(size(state));
    rho = state(1:3);
    rho_dot = state(4:6);
    r0 = state(7);
    r0_dot = state(8);
    theta0 = state(9);
    theta0_dot = state(10);

    statedot(1:3) = rho_dot;
    statedot(7) = r0_dot;
    statedot(9) = theta0_dot;

    r0_ddot = r0 * theta0_dot^2 - mu / r0^2;
    statedot(8) = r0_ddot;
    theta0_ddot = - 2 * r0_dot * theta0_dot / r0;
    statedot(10) = theta0_ddot;
    
    denom = (sqrt((r0 + rho(1))^2 + rho(2)^2 + rho(3)^2))^3;
    statedot(4) = 2 * theta0_dot * state(5) + theta0_ddot * rho(2) + theta0_dot^2 * rho(1) - mu * (r0 + rho(1)) / denom + mu / r0^2;
    statedot(5) = - 2 * theta0_dot * state(4) - theta0_ddot * rho(1) + theta0_dot^2 * rho(2) - mu * rho(2) / denom;
    statedot(6) = - mu * rho(3) / denom;
end

function statedot = NERM_maneuver(t, state, mu)
    statedot = zeros(size(state));
    rho = state(1:3);
    rho_dot = state(4:6);
    r0 = state(7);
    r0_dot = state(8);
    theta0 = state(9);
    theta0_dot = state(10);
    state_chief = state(11:16);
    state_deputy = state(17:22);

    statedot(1:3) = rho_dot;
    statedot(7) = r0_dot;
    statedot(9) = theta0_dot;
    statedot(11:13) = state(14:16);
    statedot(17:19) = state(20:22);

    statedot(14:16) = - mu * state(11:13) / norm(state(11:13))^3;
    statedot(20:22) = - mu * state(17:19) / norm(state(17:19))^3;

    r0_ddot = r0 * theta0_dot^2 - mu / r0^2;
    statedot(8) = r0_ddot;
    theta0_ddot = - 2 * r0_dot * theta0_dot / r0;
    statedot(10) = theta0_ddot;
    
    denom = (sqrt((r0 + rho(1))^2 + rho(2)^2 + rho(3)^2))^3;
    statedot(4) = 2 * theta0_dot * state(5) + theta0_ddot * rho(2) + theta0_dot^2 * rho(1) - mu * (r0 + rho(1)) / denom + mu / r0^2;
    statedot(5) = - 2 * theta0_dot * state(4) - theta0_ddot * rho(1) + theta0_dot^2 * rho(2) - mu * rho(2) / denom;
    statedot(6) = - mu * rho(3) / denom;
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

function [a, e, i, omega, RAAN, nu] = Keplerian_elements(state, mu)
    pos = state(1:3);
    vel = state(4:6);
    h = cross(pos, vel);
    W = h / norm(h);
    i = atan2(sqrt(W(1)^2 + W(2)^2), W(3));
    RAAN = atan2(W(1), - W(2));
    p = norm(h)^2 / mu;
    a = 1 / (2 / norm(pos) - norm(vel)^2 / mu);
    n = sqrt(mu / a^3);
    e = sqrt(1 - p / a);
    E = atan2((dot(pos, vel) / (a^2 * n)), (1 - norm(pos)/a));
    nu = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    u = atan2((pos(3) / sin(i)), (pos(1) * cos(RAAN) + pos(2) * sin(RAAN)));
    omega = wrapTo2Pi(u - nu);
end