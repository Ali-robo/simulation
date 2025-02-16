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

u = zeros(n+1,1);

%% u bestimmen durch numerisches berechen der ersten Schritte
ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
              x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];

[time_sol,sol] = ode45(ode,[0 h], initial_conditions);

%u(1) =  c3 * (sol(1,3)- sol(1,1)) + d3 * (sol(1,2) - sol(1,4));
%u(2) =  c3 * (sol(end,3)- sol(end,1)) + d3 * (sol(end,2) - sol(end,4));

u(1) = c3 * (sys2(1,1) - sys1(1,1)) + d3 * (sys2(1,2) - sys1(1,2)); %Anfang fehlerbehaftet.
u(2) = u(1);

%%Numerische Berechnung des Systems

ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
              x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];

[time_sol,sol] = ode45(ode,time, initial_conditions);

%Speicher der Zeitmessung
timeToc = zeros(n-1,1); 

for index = linspace(2, n, n-1)
    tic
    
    inital1 = [sys1(index-1,1);sys1(index-1,2)];
    inital2 = [sys2(index-1,1);sys2(index-1,2)];
    

    [t,sol1] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u(index-1), u(index), t_sol-h, h),[t_sol, t_sol+h], inital1);

    [t,sol2] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u(index-1), -u(index), t_sol-h, h),[t_sol, t_sol+h], inital2);

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



function dxdt = dgl_f(t,x,v,c,d,m, u0, u1, t0,h)
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        +u0 + (t-t0) * (u1-u0)/(h))/m;
end


function dxdt = dgl_d(t,x,v,c,d,m,x0,x1,v0,v1,c3,d3,t0,h)
    
    dxdt = zeros(2,1);
    dxdt(1) = v;
    dxdt(2) = (-c * x - d * v ...
        + c3 * (interp(t,x0,x1,t0,h) - x)...
        + d3 * (interp(t,v0,v1,t0,h) - v))/m;

    function p = interp(t,p0,p1,t0,h)
        p = p0 + (t-t0) * (p1-p0)/h;
    end

end

%% Ideen
%{


Zeitpunkte der ode45 vorgeben - sol1 kann vorher initialisiert werden

Die anderen Übertragungsfunktionen implementieremn

%}