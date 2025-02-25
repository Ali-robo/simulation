function dxdt = dgl_d(t,x,v,c,d,m,ux,uv,c3,d3,t0,h)
    
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        + c3 * (interp(t,ux,t0,h) - x)...
        + d3 * (interp(t,uv,t0,h) - v))/m;

    function p = interp(t,u,t0,h)
        p = u(1) + (t-t0) * (u(2)-u(1))/h;
    end

end