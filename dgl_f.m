function dxdt = dgl_f(t,x,v,c,d,m,u,t0,k)

    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        +interpolation(t,u,t0,k))/m;
end
