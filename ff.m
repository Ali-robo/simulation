 function data = ff(n,h,sysPar,init,k)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    sys1 = zeros(n+1,2);
    sys2 = zeros(n+1,2);

    sys1(1,:) = init([1 2]);
    sys2(1,:) = init([3 4]);


    %% ersten u bestimmen
    if (k < 0) && (k > 3)

        error("k: " + k + "nicht 1, 2 oder 3")
        
    elseif k == 0
        
        u = c3 * (sys2(1,1) - sys1(1,1)) + d3 * (sys2(1,2) - sys1(1,2));
    
    else

        ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
                  x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];
    
        options =  odeset(RelTol=1e-10,AbsTol=1e-12);
        [~,sol] = ode45(ode,linspace(0,-h*(k),k+1),init,options);
    
        u(1) = c3 * (sol(end,3) - sol(end,1)) + d3 * (sol(end,4)- sol(end,2));
        u(end) = c3 * (sol(1,3) - sol(1,1)) + d3 * (sol(1,4)- sol(1,2));

        for i = 1:(k-2)
            u(k-i) = c3 * (sol(i+1,3) - sol(i+1,1)) + d3 * (sol(i+1,4)- sol(i+1,2));
        end

    end

    t0 = zeros(k,1);

    %% sim
    for T = 1:n

        t0 = linspace(-h*k,0,k+1);

        [~,~,u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,h,t0,u,init);

        u = circshift(u,[0 -1]);

        u(end) = u_now;

        t0 = linspace(-h*(k-1),h,k+1);

        [sys1(T+1,:),sys2(T+1,:),u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,h,t0,u,init);

        u(end) = u_now;

        init = [sys1(T+1,:),sys2(T+1,:)];

    end

    data = struct("x1", sys1(:,1), "v1", sys1(:,2), "x2", sys2(:,1), "v2", sys2(:,2));

    data.("uDebug") = uDebug;
    data.("uPre") = uPreCalc;

end



function [sys1,sys2,u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,h,t0,u,inital)


    options =  odeset(RelTol=1e-10,AbsTol=1e-12);
    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u,t0,h), [0,h], inital([1 2]),options);% über [0,h] nachdenken
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u,t0,h), [0,h], inital([3 4]),options);
    sys2 = sol(end,[1 2]);

    u_now = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));
end