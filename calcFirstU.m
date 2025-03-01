function u = calcFirstU(sysPar,init,h)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
              x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];

    options =  odeset(RelTol=1e-10,AbsTol=1e-12);
    [~,sol] = ode45(ode,[0,-h], init,options);

    u(1) = c3 * (sol(end,3) - sol(end,1)) + d3 * (sol(end,4)- sol(end,2));
    u(2) = c3 * (sol(1,3) - sol(1,1)) + d3 * (sol(1,4)- sol(1,2));
    
end