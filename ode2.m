function dxdt = ode2(t, x, u, tiempo)
    u_int = interp1(tiempo, u, t);
    %dxdt = zeros(2,1);
    %dxdt(1) = -0.2*sqrt(x(1)-x(2)) + 0.001*u_int;
    %dxdt(2) = -0.1*sqrt(x(2)) + 0.2*sqrt(x(1)-x(2)); 
    dxdt = [-0.2*sqrt(x(1)-x(2)) + 0.001*u_int;
            -0.1*sqrt(x(2)) + 0.2*sqrt(x(1)-x(2))];
end

