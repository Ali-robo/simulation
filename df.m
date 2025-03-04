function data = df(n,h,sysPar,init)
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

    sys2_0 = [sol(end,3),sol(end,4)];


%% main loop
    for T = 1:n

        [~,temp] = ode45(@(t,x) dgl_d(t,x(1),x(2),c1,d1,m1,sys2_0,sys2(T,:),c3,d3,-h,h),[0 h], init([1 2]),options);
        sys1(T+1,:) = temp(end,:);
        [~,temp] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u(1),-u(2),-h,h),[0 h], init([3 4]));
        sys2(T+1,:) = temp(end,:);

        sys2_0 = sys2(T,:);
        u = [u(2),calcU(c3,d3,sys1(T,:),sys2(T,:))];

        [~,temp] = ode45(@(t,x) dgl_d(t,x(1),x(2),c1,d1,m1,sys2_0,sys2(T+1,:),c3,d3,0,h),[0 h], init([1 2]),options);
        sys1(T+1,:) = temp(end,:);
        [~,temp] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u(1),-u(2),0,h),[0 h], init([3 4]));
        sys2(T+1,:) = temp(end,:);

        u(2) = calcU(c3,d3,sys1(T+1,:),sys2(T+1,:));

        init = [sys1(T+1,:), sys2(T+1,:)];

    end


    data = struct("x1", sys1(:,1), "v1" , sys1(:,2), "x2", sys2(:,1), "v2", sys2(:,2));

end


function u = calcU(c3,d3,sys1,sys2)

    u = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));

end