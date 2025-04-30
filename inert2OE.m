function [a, e, i, omega, Omega, M] = inert2OE(state, mu)
r_vec = state(1:3);
v_vec = state(4:6);

v = norm(v_vec);
r = norm(r_vec);

h_vec = cross(r_vec, v_vec);
h = norm(h_vec);
nVec = cross([0, 0, 1], h_vec);
n = norm(nVec);
e_vec = (1/mu)*((v^2 - mu/r)*r_vec - dot(r_vec, v_vec)*v_vec);
e = norm(e_vec);

mech_energy = 0.5*v^2 - mu/r;
if e ~= 1
    a = -mu/(2* mech_energy);
else
    a = inf;
end

i = acos(h_vec(3)/h);
Omega = acos(nVec(1)/n);
omega = acos(dot(nVec, e_vec)/(n*e));
nu = real(acos(dot(e_vec, r_vec)/(e*r)));

if nVec(2) < 0
    Omega = 2*pi - Omega;
end
if e_vec(3) < 0
    omega = 2*pi - omega;
end
if dot (r_vec, v_vec) < 0
    nu = 2*pi - nu;
end

if i == 0 && e ~= 0
    ang = acos(e_vec(1)/e);
    if e_vec(2) < 0
        ang = 2*pi - ang;
    end
elseif i ~= 0 && e == 0
    ang = acos(dot(nVec, r_vec)/(n*r));
    if r_vec(3) < 0
        ang = 2*pi - ang;
    end
elseif i == 0 && e == 0
    ang = acos(r_vec(1)/r);
    if r_vec(2) < 0
        ang = 2*pi - ang;
    end
else
    ang = NaN;
end

M = true2mean(nu, e);
% 
% oe.a = a;
% oe.e = e;
% oe.i = i;
% oe.Omega = Omega;
% oe.omega = omega;
% oe.nu = nu;
% oe.ang = ang;
end
