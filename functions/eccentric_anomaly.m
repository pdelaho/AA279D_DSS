% eccentric_anomaly computes the eccentric anomaly using the Newton-Raphson
% method, knowing the mean anomaly and eccentricity of the spacecraft given 
% a tolerance epsilon.
%
%   Inputs:
%       M - mean anomaly [rad]
%       e - eccentricity [-]
%       epsilon - tolerance [-]
%
%   Outputs:
%       E - eccentric anomaly [rad]

function E = eccentric_anomaly(M, e, epsilon)

    M = wrapTo2Pi(M);
    E = pi;
    while abs(- E + e * sin(E) + M) / (1 - e * cos(E)) > epsilon
        E_new = E - (E - e * sin(E) - M) / (1 - e * cos(E));
        E = E_new;
    end

end

