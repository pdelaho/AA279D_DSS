% OETOECI Converts a classical set of Keplerian elements to the position and
% velocity in the Earth Centered Inertial frame
%
%   Inputs:
%       a - semi-major axis [km]
%       e - eccentricity [-]
%       inc - inclination [rad]
%       omega - argument of periapsis [rad]
%       RAAN - right ascension of the ascending node [rad]
%       nu - true anomaly [rad]
%       mu - standard gravitational parameter [km^3/s^2]
%
%   Outputs:
%       position_ECI - (3x1) position vector in ECI [km]
%       velocity_ECI - (3x1) velocity vector in ECI [km/s]

function [position_ECI, velocity_ECI] = OEtoECI(a, e, inc, omega, RAAN, nu, mu)

    % Useful quantities
    cosE = (e + cos(nu)) / (1 + e * cos(nu));
    sinE = sin(nu) * sqrt(1 - e^2) / (1 + e * cos(nu));
    n = sqrt(mu / a^3);
    
    % Position and velocity vectors in the perifocal frame
    position_perifocal = [a * (cosE - e), a * sqrt(1 - e^2) * sinE, 0];
    velocity_perifocal = a * n / (1 - e * cosE) * [-sinE, sqrt(1 - e^2) * cosE, 0];
    
    % Rotation matrix between perifocal and inertial frame
    rotation_RAAN = [cos(RAAN) sin(-RAAN) 0;
               -sin(-RAAN) cos(RAAN) 0;
               0 0 1];
    rotation_i = [1 0 0;
            0 cos(inc) sin(-inc);
            0 -sin(-inc) cos(inc)];
    rotation_omega = [cos(omega) sin(-omega) 0;
                -sin(-omega) cos(omega) 0;
                0 0 1];
    rotation_perifocalTOinertial = rotation_RAAN * rotation_i * rotation_omega;
    
    % Rotating the position and velocity vectors
    position_ECI = rotation_perifocalTOinertial * position_perifocal';
    velocity_ECI = rotation_perifocalTOinertial * velocity_perifocal';
    
end