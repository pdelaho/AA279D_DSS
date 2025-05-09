%% Newton-Raphson method to determine the optimal true anomaly for a delta e dominant case
mu = 398600.435436; % km^3/s^2

f(0.8967, 0.5, sqrt(mu / 15000^3), 307.646, 520.977*0.5)
nu = NR(pi-acos(0.5), 0.5, sqrt(mu / 15000^3), 307.646, 520.977*0.5)

function x = NR(init_guess, e, n, delta_ex_des, delta_ey_des)
%     x_prev = inf;
    x = init_guess;
    i = 0;
    f_val = f(x, e, n, delta_ex_des, delta_ey_des);

    while abs(f_val) > 1e-12 && i < 100
%         f(x, e, n, delta_ex_des, delta_ey_des)
%         f_der(x, e, n, delta_ex_des, delta_ey_des)
        new_x = x - f(x, e, n, delta_ex_des, delta_ey_des) / f_der(x, e, n, delta_ex_des, delta_ey_des);
%         x_prev = x;
        x = new_x;
        f_val = f(x, e, n, delta_ex_des, delta_ey_des);
        i = i+1;
    end
%     i
end

%% Functions

function x = f(nu, e, n, delta_ex_des, delta_ey_des)
    delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des);
    delta_ex_tild = sqrt(1-e^2) / n * (sin(nu) * sqrt(1 - delta_v^2) + ((2+e*cos(nu)) * cos(nu) + e) / (1+e*cos(nu)) * delta_v);
    delta_ey_tild = sqrt(1-e^2) / n * (-cos(nu) * sqrt(1 - delta_v^2) + ((2+e*cos(nu)) * sin(nu)) / (1+e*cos(nu)) * delta_v);
    x = delta_ey_tild - delta_ey_des / delta_ex_des * delta_ex_tild;
end

function delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des)
    if sign(delta_ex_des) == sign(delta_ey_des)
        delta_v = - sqrt(0.5 - f1(nu, e) / (2 * sqrt(4 + f1(nu, e)^2)));
    else
        delta_v = sqrt(0.5 + f1(nu, e) / (2 * sqrt(4 + f1(nu, e)^2)));
    end
end

function num = f1(nu, e)
    num = f2(nu, e) / ((1  + e * cos(nu)) * e * sin(nu));
end

function x = f2(nu, e)
    x = 2 * e^2 * cos(nu)^2 + 6 * e * cos(nu) + e^2 + 3;
end

function x = f_der(nu, e, n, delta_ex_des, delta_ey_des)
    der_delta_ex_tild_num = der_delta_ex_tild(nu, e, n, delta_ex_des, delta_ey_des);
    der_delta_ey_tild_num = der_delta_ey_tild(nu, e, n, delta_ex_des, delta_ey_des);
    x = der_delta_ey_tild_num - delta_ey_des / delta_ex_des * der_delta_ex_tild_num;
end

function x = der_delta_ex_tild(nu, e, n, delta_ex_des, delta_ey_des)
    delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des);
    der_delta_v = der_delta_vt(nu, e, delta_ex_des, delta_ey_des);
    x = sqrt(1-e^2) / n * (cos(nu) * sqrt(1-delta_v^2) - sin(nu) * der_delta_v * delta_v * (1-delta_v^2) ...
        + (delta_v * (-sin(nu) * (2+e*cos(nu)) - cos(nu) * sin(nu) * e) + der_delta_v * ((2+e*cos(nu)) * cos(nu)+e))*(1+e*cos(nu)) / (1+e*cos(nu))^2 ...
        + delta_v * ((2+e*cos(nu))*cos(nu)+e) * e * sin(nu) / (1+e*cos(nu))^2);
end

function x = der_delta_ey_tild(nu, e, n, delta_ex_des, delta_ey_des)
    delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des);
    der_delta_v = der_delta_vt(nu, e, delta_ex_des, delta_ey_des);
    x = sqrt(1-e^2) / n * (sin(nu)*sqrt(1-delta_v^2) + cos(nu)*der_delta_v*delta_v*(1-delta_v^2) ...
        + (der_delta_v*sin(nu)*(2+e*cos(nu)) +  delta_v*(cos(nu) * (2+e*cos(nu)) - sin(nu)^2*e)) * (1+e*cos(nu)) / (1+e*cos(nu))^2 ...
        + (2+e*cos(nu))*sin(nu)*delta_v*e*sin(nu) / (1+e*cos(nu))^2);
end

function x = der_delta_vt(nu, e, delta_ex_des, delta_ey_des)
    delta_v = delta_vt(nu, e, delta_ex_des, delta_ey_des);
    f1_num = f1(nu, e);
    der_f1_num = der_f1(nu, e);
    x = 0.5 * (-1 / delta_v) * (der_f1_num * 2 * (4+f1_num^2)^(0.5) / (4 * (4 +f1_num^2)) ...
        - f1_num^2 * der_f1_num * (4+f1_num^2)^(-0.5) / (4*(4+f1_num^2)));
end

function x = der_f1(nu, e)
    f2_num = f2(nu, e);
    der_f2_num = der_f2(nu, e);
    x = (der_f2_num * (1+e*cos(nu)) * e * sin(nu) - f2_num * (e*cos(nu)*(1+e*cos(nu)) - e^2*sin(nu)^2)) / ((1+e*cos(nu))^2 *e^2 * sin(nu)^2);
end

function x = der_f2(nu, e)
    x = -4*e^2*cos(nu)*sin(nu) - 6*e*sin(nu);
end