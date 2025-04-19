% E2F Computes the true anomaly from the eccentric anomaly
%   Inputs:
%       E - eccentric anomaly [rad]
%       e - eccentricity
%
%   Output:
%       f - true anomaly [rad]

function f = E2f(E, e)

    cos_f = (cos(E) - e) / (1 - e * cos(E));
    sin_f = sqrt(1 - e^2) * sin(E) / (1 - e * cos(E));
    f = atan2(sin_f, cos_f);

end