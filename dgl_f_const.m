function dxdt = dgl_f_const(t,x,v,c,d,m,u)

    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v + u * t)/m;

end