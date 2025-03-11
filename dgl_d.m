function dxdt = dgl_d(t,x,v,c,d,m,sys0,c3,d3,t0,k)
    
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        + c3 * (interpolation(t,sys0(:,1)',t0,k) - x)...
        + d3 * (interpolation(t,sys0(:,2)',t0,k) - v))/m;
end
