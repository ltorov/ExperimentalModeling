function dxdt = ode3(t, x, u, tiempo)
    u_int = interp1(tiempo, u, t);
    dxdt = [-0.1*x(1)*x(2) + u_int; 0.1*x(1)*x(2) - 0.9*x(2)];
end

