function data = ff(n,h,sysPar,init)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    
    u = calcFirstU(sysPar,init,h);

    sys1 = zeros(n+1,2);
    sys2 = zeros(n+1,2);


    sys1(1,:) = init([1 2]);
    sys2(1,:) = init([3 4]);

    uDebug = zeros(n+2,1);
    uDebug([1 2]) = u;

    for index = 1:n

        [sys1(index+1,:),sys2(index+1,:),u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,-h,h,u,init); %%vorher 0

        u = [u(2),u_now];

        [sys1(index+1,:),sys2(index+1,:),u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,0,h,u,init);

        u(2) = u_now;

        uDebug(index + 2) = u_now;

        init = [sys1(index+1,:),sys2(index+1,:)];

    end

    data = struct("x1", sys1(:,1), "v1", sys1(:,2), "x2", sys2(:,1), "v2", sys2(:,2));
    
    figure
    size(uDebug)
    plot(linspace(-h,h*n,n+2),uDebug); hold on;

end



function [sys1,sys2,u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t0,h,u, inital)

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u(1),u(2),t0,h), [0,h], inital([1 2]));% über [0,h] nachdenken
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u(1),-u(2),t0,h), [0,h], inital([3 4]));
    sys2 = sol(end,[1 2]);

    u_now = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));
end