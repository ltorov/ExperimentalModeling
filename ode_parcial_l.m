function dxdt = ode_parcial_l(t, x, u1, u2, tiempo, odes)
    odes = matlabFunction(odes);    
    u_int1 = interp1(tiempo, u1, t);
    u_int2 = interp1(tiempo, u2, t);
    dxdt = odes(u_int1,u_int2,x(1),x(2));
end