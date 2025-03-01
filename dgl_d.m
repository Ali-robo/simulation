function dxdt = dgl_d(t,x,v,c,d,m,sys0,sys1,c3,d3,t0,h)
    
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        + c3 * (interp(t,[sys0(1),sys1(1)],t0,h) - x)...
        + d3 * (interp(t,[sys0(2),sys1(2)],t0,h) - v))/m;

    function p = interp(t,u,t0,h)
        p = u(1) + (t-t0) * (u(2)-u(1))/h;
    end

end