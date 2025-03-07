 function data = ff(n,h,sysPar,init,k)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    sys1 = zeros(n+1,2);
    sys2 = zeros(n+1,2);

    sys1(1,:) = init([1 2]);
    sys2(1,:) = init([3 4]);

    options =  odeset(RelTol=1e-10,AbsTol=1e-12);


    %% ersten u bestimmen

    [u,~,~] = calcFirstU(sysPar,init,h,k);

    time = 0;
    tic;
    %% sim
    for T = 1:n

        if(mod(T,n/100) == 0)
            time = time + toc;
            disp("Noch " + ((time/(T*60))*(n-T))+ " Min")
            tic;
        end

        t0 = linspace(-h*k,0,k+1);

        [~,~,u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,k,t0,u,h,init);

        u = circshift(u,[0 -1]);

        u(end) = u_now;

        t0 = linspace(-h*(k-1),h,k+1);

        [sys1(T+1,:),sys2(T+1,:),u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,k,t0,u,h,init);

        u(end) = u_now;

        init = [sys1(T+1,:),sys2(T+1,:)];

    end

    disp("Benötigte Zeit: " + time);

    data = struct("x1", sys1(:,1), "v1", sys1(:,2), "x2", sys2(:,1), "v2", sys2(:,2));

end



function [sys1,sys2,u_now] = ode45_F_F(c1,c2,c3,d1,d2,d3,m1,m2,k,t0,u,h,inital)


    options =  odeset(RelTol=1e-10,AbsTol=1e-12);
    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u,t0,k), [0,h], inital([1 2]),options);
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u,t0,k), [0,h], inital([3 4]),options);
    sys2 = sol(end,[1 2]);

    u_now = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));
end