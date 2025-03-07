function dxdt = dgl_f(t,x,v,c,d,m,u,t0,k)

    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        +interpolation(t,u,t0,k))/m;
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




function u = lagranePoly(uIn,t0,k)

        u = @(t) 0;

        for i = 1:k

            for j = 1:k

                l = @(t) 0;

                if j ~= k

                    l = @(t) l(t) * (t - t0(j))/(t0(i)-t0(j));
                    
                end

            end

            u = @(t) u(t) + uIn(i) * l(t);
        end

end
