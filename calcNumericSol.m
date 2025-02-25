function data = calcNumericSol(sysPar,init,h,n)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    time = linspace(0,h*n,n+1);

    ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
              x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];

    [~,sol] = ode45(ode,time, init);

    data = struct("time", time,"x1", sol(:,1),"v1", sol(:,2), "x2" , sol(:,3), "v2" , sol(:,4));

end