function dxdt = dgl_d(t,x,v,c,d,m,sys0,c3,d3,t0,k)
    
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        + c3 * (interpolation(t,sys0(:,1),t0,k) - x)...
        + d3 * (interpolation(t,sys0(:,2),t0,k) - v))/m;
end



function u = interpolation(t,uIn,t0,k)

    if k == 0

        u = uIn;

    elseif k == 1

        u = uIn(1) + (t-t0(1)) * (uIn(2)-uIn(1))/(t0(2)-t0(1));

    %elseif k == 2


    else

        const = polyfit(t0,uIn,k);
        u = polyval(const,t);

    end

end