clear;
%Konstanten

par = struct( ...
    'c1', 400, ...
    'c2', 300, ...
    'c3', 80, ...
    'd1', 0.2, ...
    'd2', 0.3, ...
    'd3', 15, ...
    'm1', 6, ...
    'm2',6, ...
    'x1_0', 4, ...
    'x2_0', 2, ...
    'v1_0', 1, ...
    'v2_0',3 ...
    );


c1 = 400; d1 = 0.2; m1 = 6; c2 = 300; d2 = 0.3; m2 = 6; c3 = 80; d3 = 15;

% Anfangsbedingungen: [x1(0), v1(0), x2(0), v2(0)]
initial_conditions = [4, 2, 1, 3];


n = 100000; %Anzahl Zeitschritte am Ende
h = 0.001; %Zeitschritte
t_sol = 0; %Zeipunkt, der in der Schleife berechnent wird


%Zeitschritte für den Plot und die "richtige" Lösung
time = linspace(0,h*(n-1),n);


%System 1 und 2: [x,v] mit x,v = Spaltenvektor
sys1 = zeros(n,2);
sys2 = zeros(n,2);

%Anfangsbedingung initialisieren

sys1(1,:) = initial_conditions([1,2]);
sys2(1,:) = initial_conditions([3,4]);

%{
sys1 = struct( ...
    'x', zeros(n,1), ...
    'v', zeros(n,1));

sys2 = struct( ...
    'x', zeros(n,1), ...
    'v', zeros(n,1));

sys1.x(1) = par.x1_0; sys1.v(1) = par.v1_0;
sys2.x(1) = par.x2_0; sys2.v(1) = par.v2_0;
%}

u = zeros(n+1,1);
u(1) = c3 * (sys2(1,1) - sys1(1,1)) + d3 * (sys2(1,2) - sys1(1,2));
u(2) = u(1);

%%Numerische Berechnung des Systems

ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
              x(3); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];

[time_sol,sol] = ode45(ode,time, initial_conditions);


%%analytische Berechnung des Systems
%{
syms X1(t)  V1(t) X2(t) V2(t);

dgl_a = [diff(X1,1,t) == V1; diff(V1,1,t) == (-c1 * X1 -d1 * V1 + c3*(X2-X1) + d3 * (V2 - V1))/m1;
        diff(V2,1,t) == V2; diff(V2,1,t) == (-c2 * X2 -d2 * V2 - c3*(X2-X1) - d3 * (V2 - V1))/m2];

cond = [X1(0) == par.x1_0; V1(0) == par.v1_0; X2(0) == par.x2_0; V2(0) == par.v2_0];

sol_a = dsolve(dgl_a, cond);
%}

for index = linspace(2, n, n-1)

    inital1 = [sys1(index-1,1);sys1(index-1,2)];
    inital2 = [sys2(index-1,1);sys2(index-1,2)];
    
    dgl1 = @(t,x,v) [v; (-d1 * v - c1 * x + u(index-1) + (t-t_sol) * (u(index) - u(index-1))/h)/m1];
    dgl2 = @(t,x,v) [v; (-d2 * v - c2 * x - u(index-1) - (t-t_sol) * (u(index) - u(index-1))/h)/m2];

    [t,sol1] = ode45(@(t,x) dgl1(t,x(1),x(2)),[t_sol, t_sol+h], inital1);
    [t,sol2] = ode45(@(t,x) dgl2(t,x(1),x(2)),[t_sol, t_sol+h], inital2);

    u(index+1) = c3 * (sol2(end,1) - sol1(end,1)) + d3 * (sol2(end,2)- sol1(end,2));

    dgl1 = @(t,x,v) [v; (-d1 * v - c1 * x + u(index) + (t-t_sol) * (u(index+1) - u(index))/h)/m1];
    dgl2 = @(t,x,v) [v; (-d2 * v - c2 * x - u(index) - (t-t_sol) * (u(index+1) - u(index))/h)/m2];

    [t,sol1] = ode45(@(t,x) dgl1(t,x(1),x(2)),[t_sol, t_sol+h], inital1);
    [t,sol2] = ode45(@(t,x) dgl2(t,x(1),x(2)),[t_sol, t_sol+h], inital2);
    
    sys1(index, 1) = sol1(end,1); sys1(index, 2) = sol1(end,2);
    sys2(index, 1) = sol2(end,1); sys2(index, 2) = sol2(end,2);

    u(index+1) = c3 * (sys2(index, 1) - sys1(index, 1)) + d3 * (sys2(index, 2) - sys1(index, 2));

    t_sol = t_sol + h;

     if(mod(index,(n/100))==0)
        fprintf("|");
    end
end



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



%% Ideen


%{
Zeitpunkte der ode45 vorgeben - sol1 kann vorher initialisiert werden

functions verwenden um u und dgl aufzustellen


%}