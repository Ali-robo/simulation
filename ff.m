 function data = ff(n,h,sysPar,init)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    sys1 = zeros(n+1,2);
    sys2 = zeros(n+1,2);

    sys1(1,:) = init([1 2]);
    sys2(1,:) = init([3 4]);

    u = zeros(2,1);

    %% ersten u bestimmen
    ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
              x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];

    options =  odeset(RelTol=1e-10,AbsTol=1e-12);
    [~,sol] = ode45(ode,[0,-h], init,options);

    u(1) = c3 * (sol(end,3) - sol(end,1)) + d3 * (sol(end,4)- sol(end,2));
    u(2) = c3 * (sol(1,3) - sol(1,1)) + d3 * (sol(1,4)- sol(1,2));

    uDebug = zeros(n+2,1);
    uDebug([1 2]) = u;
     
    uPreCalc = zeros(n,1);

    %% sim
    for t = 1:n

        [~,~,u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,-h,h,u,init);

        u = [u(2),u_now];
        uPreCalc(t) = u_now;

        [sys1(t+1,:),sys2(t+1,:),u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,0,h,u,init);

        u(2) = u_now;

        uDebug(t + 2) = u_now;

        init = [sys1(t+1,:),sys2(t+1,:)];

    end

    data = struct("x1", sys1(:,1), "v1", sys1(:,2), "x2", sys2(:,1), "v2", sys2(:,2));

    data.("uDebug") = uDebug;
    data.("uPre") = uPreCalc;

end



function [sys1,sys2,u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t0,h,u, inital)


    options =  odeset(RelTol=1e-10,AbsTol=1e-12);
    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u(1),u(2),t0,h), [0,h], inital([1 2]),options);% über [0,h] nachdenken
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u(1),-u(2),t0,h), [0,h], inital([3 4]),options);
    sys2 = sol(end,[1 2]);

    u_now = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));
end