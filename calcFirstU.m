function [u,sys1,sys2] = calcFirstU(sysPar,init,h,k)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    u = zeros(1,k+1);
    sys1 = zeros(k+1,2);
    sys2 = zeros(k+1,2);

    if (k < 0)

        error("k muss größer gleich 0 sein")
        
    elseif k == 0
        
        u = c3 * (init(3) - init(1)) + d3 * (init(4) - init(2));
        sys1 = init([1 2]);
        sys2 = init([3 4]);


    else

        ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
                  x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];
    
        options =  odeset(RelTol=1e-10,AbsTol=1e-12);
        [~,sol] = ode45(ode,linspace(0,-h*(k),k+1),init,options);
    
        u(1) = c3 * (sol(end,3) - sol(end,1)) + d3 * (sol(end,4)- sol(end,2));
        u(end) = c3 * (sol(1,3) - sol(1,1)) + d3 * (sol(1,4)- sol(1,2));

        sys1(1,:) = sol(end,[1 2]);
        sys1(end,:) = sol(1,[1 2]);

        sys2(1,:) = sol(end,[3 4]);
        sys2(end,:) = sol(1,[3 4]);


        for i = 1:(k-2)
            u(k-i) = c3 * (sol(i+1,3) - sol(i+1,1)) + d3 * (sol(i+1,4)- sol(i+1,2));

            sys1(k-i,:) = sol(i+1,[1 2]);
            sys2(k-i,:) = sol(i+1,[3 4]);

        end

    end
end