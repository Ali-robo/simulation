clear;

%Konstanten

par = struct( ...
    'c1', 400, ...
    'c2', 300, ...
    'cc', 80, ...
    'd1', 0.2, ...
    'd2', 0.3, ...
    'dc', 15, ...
    'm1', 6, ...
    'm2',6, ...
    'x1_0', 4, ...
    'x2_0', 2, ...
    'v1_0', 1, ...
    'v2_0',3 ...
    );



n = 100000; %Anzahl Zeitschritte am Ende
h = 0.001; %Zeitschritte
t_sol = 0; %Zeipunkt, der in der Schleife berechnent wird


%Zeitschritte für den Plot und die "richtige" Lösung
time = linspace(0,h*(n-1),n);


% Anfangsbedingungen: [x1(0), v1(0), x2(0), v2(0)]
initial_conditions = [par.x1_0, par.v1_0, par.x2_0, par.v2_0];


sys1 = struct( ...
    'x', zeros(n,1), ...
    'v', zeros(n,1));

sys2 = struct( ...
    'x', zeros(n,1), ...
    'v', zeros(n,1));

sys1.x(1) = par.x1_0; sys1.v(1) = par.v1_0;
sys2.x(1) = par.x2_0; sys2.v(1) = par.v2_0;

u = zeros(n+1,1); u(1) = par.cc * (par.x2_0 - par.x1_0) + par.dc * (par.v2_0 - par.v1_0);
                  u(2) = u(1);

%%Numerische Berechnung des Systems

ode = @(t,x) [x(2); (-par.c1 * x(1) -par.d1 * x(2) + par.cc*(x(3)-x(1)) + par.dc * (x(4) - x(2)))/par.m1;
              x(3); (-par.c2 * x(3) -par.d2 * x(4) - par.cc*(x(3)-x(1)) - par.dc * (x(4) - x(2)))/par.m2];

[time_sol,sol] = ode45(ode,time, initial_conditions);


%%analytische Berechnung des Systems
%{
syms X1(t)  V1(t) X2(t) V2(t);

dgl_a = [diff(X1,1,t) == V1; diff(V1,1,t) == (-par.c1 * X1 -par.d1 * V1 + par.cc*(X2-X1) + par.dc * (V2 - V1))/par.m1;
        diff(V2,1,t) == V2; diff(V2,1,t) == (-par.c2 * X2 -par.d2 * V2 - par.cc*(X2-X1) - par.dc * (V2 - V1))/par.m2];

cond = [X1(0) == par.x1_0; V1(0) == par.v1_0; X2(0) == par.x2_0; V2(0) == par.v2_0];

sol_a = dsolve(dgl_a, cond);
%}

for index = linspace(2, n, n-1)

    inital1 = [sys1.x(index-1);sys1.v(index-1)];
    inital2 = [sys2.x(index-1);sys2.v(index-1)];
    
    dgl1 = @(t,x,v) [v; (-par.d1 * v - par.c1 * x + u(index-1) + (t-t_sol) * (u(index) - u(index-1))/h)/par.m1];
    dgl2 = @(t,x,v) [v; (-par.d2 * v - par.c2 * x - u(index-1) - (t-t_sol) * (u(index) - u(index-1))/h)/par.m2];

    
    [t,sol1] = ode45(@(t,x) dgl1(t,x(1),x(2)),[t_sol, t_sol+h], inital1);
    [t,sol2] = ode45(@(t,x) dgl2(t,x(1),x(2)),[t_sol, t_sol+h], inital2);

    u(index+1) = par.cc * (sol2(end,1) - sol1(end,1)) + par.dc * (sol2(end,2)- sol1(end,2));

    dgl1 = @(t,x,v) [v; (-par.d1 * v - par.c1 * x + u(index) + (t-t_sol) * (u(index+1) - u(index))/h)/par.m1];
    dgl2 = @(t,x,v) [v; (-par.d2 * v - par.c2 * x - u(index) - (t-t_sol) * (u(index+1) - u(index))/h)/par.m2];

    [t,sol1] = ode45(@(t,x) dgl1(t,x(1),x(2)),[t_sol, t_sol+h], inital1);
    [t,sol2] = ode45(@(t,x) dgl2(t,x(1),x(2)),[t_sol, t_sol+h], inital2);
    
    sys1.x(index) = sol1(end,1); sys1.v(index) = sol1(end,2);
    sys2.x(index) = sol2(end,1); sys2.v(index) = sol2(end,2);

    u(index+1) = par.cc * (sys2.x(index) - sys1.x(index)) + par.dc * (sys2.v(index) - sys1.v(index));



    t_sol = t_sol + h;

     if(mod(index,(n/100))==0)
        fprintf("|");
    end
end



%% plotting

nexttile

title("System 1");

plot(time,sys1.x, "Color","r"); hold on;
plot(time,sys1.v, "Color","b");

plot(time_sol, sol(:,1), "Color","r","LineStyle","--");
plot(time_sol, sol(:,2), "Color","b","LineStyle","--"); hold off;

ylim([-100,100]);
xlabel("time");
ylabel("x1, v1");
legend("x1", "v1", "x1x", "v1n");


nexttile

title("System 2");

plot(time,sys2.x, "Color","r"); hold on;
plot(time,sys2.v, "Color","b"); 

plot(time_sol, sol(:,3), "Color","g","LineStyle","--");
plot(time_sol, sol(:,4), "Color","g","LineStyle","--"); hold off;

ylim([-100,100]);
xlabel("time");
ylabel("x2, v2");
legend("x2", "v2","x1x", "v1n");
