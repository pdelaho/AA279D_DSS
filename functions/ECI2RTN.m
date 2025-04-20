% ECItoRTN Computes the rotation matrix from the Earth Centered Inertial
% frame to the RTN frame
%   Input:
%       state = [position_ECI, velocity_ECI]
%           position_ECI - position in inertial frame
%           velocity_ECI - velocity in inertial frame
%
%   Output:
%       rotation - (3x3) rotation matrix from inertial to RTN frames

function state_RTN = ECI2RTN(state_frame, state_ECI, mu)

    % Reformatting the input into row vectors
    position = zeros(1,3);
    velocity = zeros(1,3);
    position(:) = state_frame(1:3);
    velocity(:) = state_frame(4:6);
    state_ECI_row = reshape(state_ECI, 1, []);

    % RTN unit vectors
    R = position / norm(position);
    h = cross(position, velocity);
    N = h / norm(h);
    T = cross(N, R);

    % Rotation matrix from ECI to RTN
    rotation = [R; T; N];

    [a, e, i, omega, RAAN, f] = ECI2OE_f(state_frame, mu);

    theta_dot_vec = [0 0 sqrt(mu / (a^3 * (1 - e^2)^3)) * (1 + e * cos(f))^2];

    state_RTN = zeros(1,6);
    state_RTN(1:3) = rotation * state_ECI_row(1:3)';
    state_RTN(4:6) = rotation * state_ECI_row(4:6)' - cross(theta_dot_vec, state_RTN(1:3)')';
    
end