function dxdt = ode(t, x, u, tiempo, odes)
    u_int = interp1(tiempo, u, t);
    dxdt = [double(subs(odes(1),{'x1','x2','u'},{x(1),x(2),u_int}));
            double(subs(odes(2),{'x1','x2','u'},{x(1),x(2),u_int}))];
end

