% ECI2OE_f Computes the classical Keplerian orbital elements (with the true
% anomaly) from the position and velocity of the spacecraft in the Earth 
% Centered Inertial frame.
%
%   Inputs:
%       state - position [km] and velocity [km/s] in ECI
%       mu - standard gravitational parameter [km^3/s^2]
%
%   Outputs:
%       a - semi-major axis [km]
%       e - eccentricity [-]
%       i - inclination [rad]
%       omega - argument of periapsis [rad]
%       RAAN - right ascension of the ascending node [rad]
%       f - true anomaly [rad]

function [a, e, i, omega, RAAN, f] = ECI2OE_f(state, mu)

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
    f = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2));
    u = atan2((pos(3) / sin(i)), (pos(1) * cos(RAAN) + pos(2) * sin(RAAN)));
    omega = wrapTo2Pi(u - f);
    
end