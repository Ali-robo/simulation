function dxdt = dgl_f(t,x,v,c,d,m,u,t0,h)
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        +u0 + (t-t0) * (u1-u0)/(h))/m;
end


function u = interpolation(t,uIn,t0,k)      %%polyval und polyfit

    if k == 0

        u = uIn * t;

    elseif k == 1

        u = uIn(1) + (t-t0(1)) * (uIn(2)-uIn(1))/(t0(2)-t0(1));

    elseif k == 2

        u = uIn

    elseif k == 3


    end

end