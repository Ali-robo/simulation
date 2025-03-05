function data = ff_const(n,h,sysPar,init)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;


    options =  odeset(RelTol=1e-10,AbsTol=1e-12);

    sys1 = zeros(n+1,2);
    sys2 = zeros(n+1,2);
    sys1(1,:) = init([1 2]);
    sys2(1,:) = init([3 4]);

    u = c3 * (sys2(1,1) - sys1(1,1)) + d3 * (sys2(1,2) - sys1(1,2));


    for T = 1:n


        [~,temp] = ode45(@(t,x) dgl_f_const(t,x(1),x(2),c1,d1,m1,u),[0 h], sys1(T,:),options);
        sys1(T+1,:) = temp(end,:);
        [~,temp] = ode45(@(t,x) dgl_f_const(t,x(1),x(2),c2,d2,m2,-u),[0 h], sys2(T,:),options);
        sys2(T+1,:) = temp(end,:);

        u = c3 * (sys2(T+1,1) - sys1(T+1,1)) + d3 * (sys2(T+1,2) - sys1(T+1,2));

        [~,temp] = ode45(@(t,x) dgl_f_const(t,x(1),x(2),c1,d1,m1,u),[0 h], sys1(T,:),options);
        sys1(T+1,:) = temp(end,:);
        [~,temp] = ode45(@(t,x) dgl_f_const(t,x(1),x(2),c2,d2,m2,-u),[0 h], sys2(T,:),options);
        sys2(T+1,:) = temp(end,:);

        u = c3 * (sys2(T+1,1) - sys1(T+1,1)) + d3 * (sys2(T+1,2) - sys1(T+1,2));


    end

    data = struct("x1", sys1(:,1), "v1", sys1(:,2), "x2", sys2(:,1), "v2", sys2(:,2));


end
