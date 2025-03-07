function data = dd(n,h,sysPar,init,k)

    
    

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    sys1 = zeros(n+1,2);
    sys2 = zeros(n+1,2);

    sys1(1,:) = init([1 2]);
    sys2(1,:) = init([3 4]);

    options =  odeset(RelTol=1e-10,AbsTol=1e-12);

    [~,sys1_0,sys2_0] = calcFirstU(sysPar,init,h,k);


    %% main loop

    for T = 1:n

        t0 = linspace(-h*k,0,k+1);

        [~,temp] = ode45(@(t,x) dgl_d(t,x(1),x(2),c1,d1,m1,sys2_0,c3,d3,t0,k),[0 h],sys1(T,:),options);
        sys1(T+1,:) = temp(end,:);
        [~,temp] = ode45(@(t,x) dgl_d(t,x(1),x(2),c2,d2,m2,sys1_0,c3,d3,t0,k),[0 h],sys2(T,:),options);
        sys2(T+1,:) = temp(end,:);


        sys1_0 = circshift(sys1_0,[-1 0]);
        sys2_0 = circshift(sys2_0,[-1 0]);

        sys1_0(end,:) = sys1(T,:);
        sys2_0(end,:) = sys2(T,:);

        t0 = linspace(-h*(k-1),h,k+1);

        [~,temp] = ode45(@(t,x) dgl_d(t,x(1),x(2),c1,d1,m1,sys2_0,c3,d3,t0,k),[0 h],sys1(T,:),options);
        sys1(T+1,:) = temp(end,:);
        [~,temp] = ode45(@(t,x) dgl_d(t,x(1),x(2),c2,d2,m2,sys1_0,c3,d3,t0,k),[0 h],sys2(T,:),options);
        sys2(T+1,:) = temp(end,:);

        sys1_0(end,:) = sys1(T,:);
        sys2_0(end,:) = sys2(T,:);

    end

    data = struct("x1", sys1(:,1), "v1", sys1(:,2), "x2", sys2(:,1), "v2", sys2(:,2));

end