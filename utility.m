function [ct, at, ht, ut] = utility(par, gr, b, theta, w_l, w_h, r, P_h)

amin= -w_l/(1+r) + 1.0e-2;
amax= b + w_l - 1.0e-2;
a_0 = bisection(@(x) (b - x + w_l)^(-par.psi) - par.beta*(1+r)*(w_l + (1+r)*x)^(-par.psi), amin, amax);

if a_0 < - gr.Abar
    a_0 = -gr.Abar;
elseif a_0 > b + w_l
    a_0 = b + w_l - eps;
end

if -gr.Abar < b-P_h
    amin = (-w_l/(1+r) + 1.0e-10)*(theta~=1) + (-w_h/(1+r) + 1.0e-10)*(theta==1);
    amax = b - P_h - 1.0e-10;
    a_1 = bisection(@(x) (b - x + - P_h)^(-par.psi) - par.beta*(1+r)*( (theta^par.eta)*(w_h + (1+r)*x)^(-par.psi) + (1-(theta^par.eta))*(w_l + (1+r)*x)^(-par.psi)), amin, amax);   % With college education

    % Savings must be above borrowing constraint, and cannot make consumption
    % negative
    if a_1 < - gr.Abar
        a_1 = -gr.Abar;
    elseif a_1 > b - P_h
        a_1 = b - P_h - eps;
    end

    u_0 = ((b - a_0 + w_l)^(1-par.psi))/(1-par.psi) + par.beta * ((w_l + (1+r)*a_0)^(1-par.psi))/(1-par.psi);
    u_1 = ((b - a_1 - P_h)^(1-par.psi))/(1-par.psi) + ...
            par.beta*( (theta^par.eta)*((w_h + (1+r)*a_1)^(1-par.psi))/(1-par.psi) + ...
                 (1-(theta^par.eta))*((w_l + (1+r)*a_1)^(1-par.psi))/(1-par.psi));

    if u_0 >= u_1
        ht = 0;
        at = a_0;
        ut = u_0;
        ct = b - a_0 + w_l;
    else
        ht = 1;
        at = a_1;
        ut = u_1;
        ct = b - a_1 - P_h;
    end
else
    ht = 0;
    at = a_0;
    ut = ((b - a_0 + w_l)^(1-par.psi))/(1-par.psi) + par.beta * ((w_l + (1+r)*a_0)^(1-par.psi))/(1-par.psi);
    ct = b - a_0 + w_l;
end

end