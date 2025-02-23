clear;

%Konstanten
c1 = 400; d1 = 0.2; m1 = 6; c2 = 300; d2 = 0.3; m2 = 6; c3 = 80; d3 = 15;

% Anfangsbedingungen: [x1(0), v1(0), x2(0), v2(0)]
initial_conditions = [4, 2, 1, 3];

n = 10000; %Anzahl Zeitschritte am Ende
h = 0.001; %Zeitschritte
t_sol = 0; %Zeipunkt, der in der Schleife berechnent wird

%Zeitschritte für den Plot und die "richtige" Lösung
time = linspace(0,h*(n-1),n);

%System 1 und 2: [x,v] mit x,v = Spaltenvektoren
sys1 = zeros(n,2);
sys2 = zeros(n,2);

%Anfangsbedingung initialisieren
sys1(1,:) = initial_conditions([1,2]);
sys2(1,:) = initial_conditions([3,4]);


%%Numerische Berechnung des Systems

ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
              x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];

[time_sol,sol] = ode45(ode,time, initial_conditions);


%%Numerische Berechnung des ersten Zeitschrittes

[~,sol_init] = ode45(ode,[0,-h], initial_conditions);   %sol_init = [x1,v1,x2,v2] mit index = 1 t = 0, index = end t = -h



%Speicher der Zeitmessung
%timeToc = zeros(n-1,1);

%% Force - Force

%{
initial_conditions = [initial_conditions, u(c3,d3,sol_init(end,[1 2]), sol_init(end,[3 4])), u(c3,d3,sol_init(1,[1 2]), sol_init(1,[3 4]))];

for index = linspace(1,n,n)

    [sys1(index,:), sys2(index,:),u_now] = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol-h,t_sol,h,initial_conditions);

    initial_conditions([5,6]) = [initial_conditions(6),u_now];

    [sys1(index,:), sys2(index,:),u_now] = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol,t_sol,h,initial_conditions);

    initial_conditions(6) = u_now;

    initial_conditions([1,2,3,4]) = [sys1(index,:),sys2(index,:)];

    t_sol = t_sol + h;

    if(mod(index,(n/100))==0)
        fprintf("|");
    end
end


%% Force - Displacement

initial_conditions = [initial_conditions, u(c3,d3,sol_init(end,[1 2]), sol_init(end,[3 4])), u(c3,d3,sol_init(1,[1 2]),0,0, sol_init(1,[3 4])), [sol_init(1,1),sol_init(end,1)], [sol_init(1,2),sol_init(end,2)]];

for index = linspace(1,n,n)

    [sys1(index,:), sys2(index,:),u_now] = cosim_F_D(c1,c2,c3,d1,d2,d3,m1,m2,u(1),u(2),sys1(index-2,1),sys1(index-2, 2),sys1(index-1,1),sys1(index-1, 2), t_sol-h,t_sol,h,initial_conditions);

    initial_conditions([5,6,7,8,9,10,11,12]) = [initial_conditions(6), u_now, 0, 0, initial_conditions(10), sys1(index,1), initial_conditions(12), sys1(index,2)];

    [sys1(index,:), sys2(index,:),u_now] =cosim_F_D(c1,c2,c3,d1,d2,d3,m1,m2,u(1),u(2),sys1(index-1,1),sys1(index-1, 2),sys1(index,1),sys1(index, 2), t_sol,t_sol,h,initial_conditions);

    %u(2) = u_now;

    initial_conditions = [sys1(index,:),sys2(index,:)];

    t_sol = t_sol + h;

    if(mod(index,(n/100))==0)
        fprintf("|");
    end
end
%}
%% Displacement Displacement

% = [x1_0, v1_0, x2_0, v2_0, u11_0, u11_1, u12_0, u12_1, u21_0, u21_1, u22_0, u22_1]

% = [ x1_-1, v1_-1, x2_-1, v2_-1;x1_0, v1_0, x2_0, v2_0]
initial_conditions = [sol_init(end,:);initial_conditions];
debug = initial_conditions;

for index = linspace(1,n,n)

   [sys1(index,:), sys2(index,:),u_now]  = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol-h,t_sol,h,initial_conditions);

   initial_conditions([5,6,7,8,9,10,11,12]) = [initial_conditions(6), u_now(1), initial_conditions(8), u_now(2), initial_conditions(10), u_now(3), initial_conditions(12), u_now(4)];

   [sys1(index,:), sys2(index,:),u_now] = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol,t_sol,h,initial_conditions);

   initial_conditions([1 2 3 4 6 8 10 12]) = [sys1(index,:),sys2(index,:), u_now];
    
end




%{
for index = linspace(2, 2*n, 2*n-1)
    tic
    
    inital1 = [sys1(index-1,1);sys1(index-1,2)];
    inital2 = [sys2(index-1,1);sys2(index-1,2)];
    

    [t,sol1] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1, u(index-1), u(index), t_sol-h, h),[t_sol, t_sol+h], inital1);
    [t,sol2] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u(index-1),-u(index), t_sol-h, h),[t_sol, t_sol+h], inital2);

    u(index+1) = c3 * (sol2(end,1) - sol1(end,1)) + d3 * (sol2(end,2)- sol1(end,2));

    [t,sol1] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u(index), u(index+1), t_sol, h),[t_sol, t_sol+h], inital1);
    [t,sol2] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u(index), -u(index+1), t_sol, h),[t_sol, t_sol+h], inital2);
    
    sys1(index, 1) = sol1(end,1); sys1(index, 2) = sol1(end,2);
    sys2(index, 1) = sol2(end,1); sys2(index, 2) = sol2(end,2);

    u(index+1) = c3 * (sys2(index, 1) - sys1(index, 1)) + d3 * (sys2(index, 2) - sys1(index, 2));

    t_sol = t_sol + h;
    
    timeToc(index-1) = toc;

     if(mod(index,(n/100))==0)
        fprintf("|");
    end
end

mean(timeToc)
%}

%% plotting

nexttile

title("System 1");

plot(time,sys1(:,1), "Color","r"); hold on;
plot(time,sys1(:,2), "Color","b");

plot(time_sol, sol(:,1), "Color","r","LineStyle","--");
plot(time_sol, sol(:,2), "Color","b","LineStyle","--"); hold off;

ylim([-100,100]);
xlabel("time");
ylabel("x1, v1");
legend("x1", "v1", "x1x", "v1n");


nexttile

title("System 2");

plot(time,sys2(:,1), "Color","r"); hold on;
plot(time,sys2(:,2), "Color","b"); 

plot(time_sol, sol(:,3), "Color","g","LineStyle","--");
plot(time_sol, sol(:,4), "Color","g","LineStyle","--"); hold off;

ylim([-100,100]);
xlabel("time");
ylabel("x2, v2");
legend("x2", "v2","x1x", "v1n");

%%  Functions and DGL


function dxdt = dgl_f(t,x,v,c,d,m, u0, u1, t0,h)
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        +u0 + (t-t0) * (u1-u0)/(h))/m;
end


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



function [sys1,sys2,u] = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t0,t_sol,h,inital,kopplung)  %sol = [[x1,v1],[x2,v2],u]  mit inital = [x1_0,v1_0,x2_0,v2_0,u0,u1]

    u0 = u_calc(c3,d3,inital(1,[1 2]), inital(1,[3 4]));   
    u1 = u_calc(c3,d3,inital(2,[1 2]), inital(1,[3 4]));
   

   [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u0,u1,t0,h), [t_sol, t_sol + h], inital([1 2]));
   sys1 = sol(end,[1 2]);

   [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u0,-u1,t0,h), [t_sol, t_sol + h], inital([3 4]));
   sys2 = sol(end,[1 2]);

   u = zeros(1,4);
   u(1) = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));

   function u = u_calc(c3,d3,sys1,sys2)
        
    u = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));

end

end

function [sys1,sys2,u] = cosim_F_D(c1,c2,c3,d1,d2,d3,m1,m2,t0,t_sol,h,inital)  %sol = [[x1,v1],[x2,v2],u] mit inital = [x1_0,v1_0,x2_0,v2_0,u0,u1,ux,uv] mit ux = [x2_0,x2_1] uv = [v2_0,v2_1]

   [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,inital(5),inital(6),t0,h), [t_sol, t_sol + h], inital([1 2]));
   sys1 = sol(end,[1 2]);

   [~,sol] = ode45(@(t,x) dgl_d(t,x(1),x(2),c2,d2,m2,inital(7),inital(8),c3,d3,t0,h), [t_sol,t_sol + h], inital([3 4]));
   sys2 = sol(end,[1 2]);

   u = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));

end

function [sys1,sys2,u] = cosim_D_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol,h,inital)  %sol = [[x1,v1],[x2,v2],u]   mit inital = [x1_0,v1_0,x2_0,v2_0,ux,uv,u0,u1] mit ux = [x1_0,x1_0] uv = [v1_0,v1_1]

    [~,sol] = ode45(@(t,x) dgl_d(t,x(1),x(2),c1,d1,m1,init(5),init(6),c3,d3,t0,h), [t_sol,t_sol + h], inital([3 4]));
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-inital(7),-inital(8),t0,h), [t_sol, t_sol + h], inital([3 4]));
    sys2 = sol(end,[1 2]);

    u = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));

end

function [sys1,sys2] = cosim_D_D(c1,c2,c3,d1,d2,d3,m1,m2,t_sol,h,inital)  %sol = [[x1,v1],[x2,v2]]  mit inital = [x1_0,v1_0,x2_0,v2_0,ux1,uv1,ux2,uv2] mit ux1 = [x1_0,x1_1] ...

    [~,sol] = ode45(@(t,x) dgl_d(t,x(1),x(2),c1,d1,m1,init(5),init(6),c3,d3,t0,h), [t_sol,t_sol + h], inital([3 4]));
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_d(t,x(1),x(2),c2,d2,m2,inital(7),inital(8),c3,d3,t0,h), [t_sol,t_sol + h], inital([3 4]));
    sys2 = sol(end,[1 2]);

end

function u = u(c3,d3,sys1,sys2)
        
    u = c3 * (sys2(1) - sys1(1)) + d3 * (sys2(2) - sys1(2));

end

%% Ideen
%{


Zeitpunkte der ode45 vorgeben - sol1 kann vorher initialisiert werden

Die anderen Übertragungsfunktionen implementieremn


Bei Force-Force müssen nur zwei u gespeichert werden, als initial
condition impelemtieren?
%}