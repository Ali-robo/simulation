clear;

%Konstanten
c1 = 400; d1 = 0.2; m1 = 6; c2 = 300; d2 = 0.3; m2 = 6; c3 = 80; d3 = 15;

% Anfangsbedingungen: [x1(0), v1(0), x2(0), v2(0)]
initial_conditions = [4, 2, 1, 3];

n = 1000; %Anzahl Zeitschritte am Ende
h = 0.01; %Zeitschritte
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


%% sim mit verschiedenen Schrittweiten - Fehlerberechnung

load('matlab.mat');
tToSim = 5;

stepSize = 1000;
nStart = 100;
nEnd = 10000;

nSim = linspace(nStart,nEnd,1+(nEnd-nStart)/(stepSize));


data = struct("maxFehler", zeros(4,size(nSim,2)));

for numOfSim = linspace(1,size(nSim,2),size(nSim,2))

    n = nSim(numOfSim);
    h = tToSim/n;
    initial_conditions = [sol_init(1,:)];
    u = [sol_init(end,:);initial_conditions];
    t_sol = 0;
    sys1 = zeros(n,2);
    sys2 = zeros(n,2);

    time = linspace(0,h*(n-1),n);
    [time_sol,sol] = ode45(ode,time, initial_conditions);

    for index = 1:n-1

       [sys1(index+1,:), sys2(index+1,:)]  = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol-h,t_sol,h,u,initial_conditions);
    
       u = [u([5 6 7 8]), sys1(index+1,:), sys2(index+1,:)];
    
       [sys1(index+1,:), sys2(index+1,:)] = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol,t_sol,h,u,initial_conditions);
    
       
       initial_conditions = [sys1(index+1,:),sys2(index+1,:)];
       u([5 6 7 8]) =  [sys1(index+1,:), sys2(index+1,:)];
    
       t_sol = t_sol + h;
    end

    data.("Sim" + numOfSim).sys1 = sys1;
    data.("Sim" + numOfSim).sys2 = sys2;

    data.("Sim" + numOfSim).fehler = abs([sys1,sys2] - sol);

    data.maxFehler(:,numOfSim) = mean(data.("Sim" + numOfSim).fehler);

    fprintf("Sim "+ numOfSim + "done." + newline);

end


%plot

timeSteps = tToSim./nSim;

loglog(timeSteps,data.maxFehler(1,:),"-or",timeSteps,data.maxFehler(2,:),"--or", timeSteps,data.maxFehler(3,:),"-ob", timeSteps,data.maxFehler(4,:),"--ob");
legend("x1","v1","x2","v2");
grid on;


%% Main Loop

load('matlab.mat')
u = [sol_init(end,:);initial_conditions];


for k = 1:2
    ax(k) = subplot(2,1,k);
end

debugText = strings(4,4);
tempData = zeros(n,4);

for index = 1:n-1

   [sys1(index+1,:), sys2(index+1,:)]  = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol-h,t_sol,h,u,initial_conditions);

   u = [u([5 6 7 8]), sys1(index+1,:), sys2(index+1,:)];

   temp = abs([sol(index+1,[1 2]) - sys1(index+1,:), sol(index+1,[3 4]) - sys2(index+1,:)]);

   debugText(1,:) = string(sol(index+1,:));
   debugText(2,:) = string([sys1(index+1,:),sys2(index+1,:)]);

   [sys1(index+1,:), sys2(index+1,:)] = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t_sol,t_sol,h,u,initial_conditions);

   
   initial_conditions = [sys1(index+1,:),sys2(index+1,:)];
   u([5 6 7 8]) =  [sys1(index+1,:), sys2(index+1,:)];

   tempData(index,:) = abs([sol(index+1,[1 2])- sys1(index+1,:), sol(index+1,[3 4]) - sys2(index+1,:)]) - temp;

   debugText(3,:) = string([sys1(index+1,:),sys2(index+1,:)]);
   debugText(4,:) = string(tempData(index,:));


   disp(debugText);
   newline;

   t_sol = t_sol+h;

   if(mod(index,(n/100))==0)
        fprintf("|");
    end
end

subplot(ax(1));
plot(time,tempData(:,1),"r"); hold on;
plot(time,tempData(:,2),"b"); hold off;
grid on;
subplot(ax(2))
plot(time,tempData(:,3),"r"); hold on;
plot(time,tempData(:,4),"b"); hold off;
grid on;


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



function [sys1,sys2] = cosim_F_F(c1,c2,c3,d1,d2,d3,m1,m2,t0,t_sol,h,kopplung, inital)  %  mit inital = [x1_0,v1_0,x2_0,v2_0] kopplung = [x1_0,v1_0,x2_0,v2_0,x1_1,v1_1,x2_1,v2_1]

    u0 = c3 * (kopplung(3)-kopplung(1)) + d3 * (kopplung(4)-kopplung(2));
    u1 = c3 * (kopplung(7)-kopplung(5)) + d3 * (kopplung(8)-kopplung(6));

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u0,u1,t0,h), [t_sol, t_sol + h], inital([1 2]));
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u0,-u1,t0,h), [t_sol, t_sol + h], inital([3 4]));
    sys2 = sol(end,[1 2]);

end


function [sys1,sys2] = cosim_F_D(c1,c2,c3,d1,d2,d3,m1,m2,t0,t_sol,h,kopplung,inital)  % mit inital = [x1_0,v1_0,x2_0,v2_0] kopplung = [x1_0,v1_0,x2_0,v2_0,x1_1,v1_1,x2_1,v2_1]

   u0 = c3 * (kopplung(3)-kopplung(1)) + d3 * (kopplung(4)-kopplung(2));
   u1 = c3 * (kopplung(7)-kopplung(5)) + d3 * (kopplung(8)-kopplung(6));

   [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c1,d1,m1,u0,u1,t0,h), [t_sol, t_sol + h], inital([1 2]));
   sys1 = sol(end,[1 2]);

   [~,sol] = ode45(@(t,x) dgl_d(t,x(1),x(2),c2,d2,m2,kopplung([1,5]),kopplung([2,6]),c3,d3,t0,h), [t_sol,t_sol + h], inital([3 4]));
   sys2 = sol(end,[1 2]);

end

function [sys1,sys2] = cosim_D_F(c1,c2,c3,d1,d2,d3,m1,m2,t0,t_sol,h,kopplung,inital)  % mit inital = [x1_0,v1_0,x2_0,v2_0] kopplung = [x1_0,v1_0,x2_0,v2_0,x1_1,v1_1,x2_1,v2_1]

    u0 = c3 * (kopplung(3)-kopplung(1)) + d3 * (kopplung(4)-kopplung(2));
    u1 = c3 * (kopplung(7)-kopplung(5)) + d3 * (kopplung(8)-kopplung(6));


    [~,sol] = ode45(@(t,x) dgl_d(t,x(1),x(2),c1,d1,m1,kopplung([3,7]),kopplung([4,8]),c3,d3,t0,h), [t_sol,t_sol + h], inital([1 2]));
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_f(t,x(1),x(2),c2,d2,m2,-u0,-u1,t0,h), [t_sol, t_sol + h], inital([3 4]));
    sys2 = sol(end,[1 2]);

end

function [sys1,sys2] = cosim_D_D(c1,c2,c3,d1,d2,d3,m1,m2,t0,t_sol,h,kopplung,inital)  % mit inital = [x1_0,v1_0,x2_0,v2_0] kopplung = [x1_0,v1_0,x2_0,v2_0,x1_1,v1_1,x2_1,v2_1]

    [~,sol] = ode45(@(t,x) dgl_d(t,x(1),x(2),c1,d1,m1,kopplung([3,7]),kopplung([4,8]),c3,d3,t0,h), [t_sol,t_sol + h], inital([1 2]));
    sys1 = sol(end,[1 2]);

    [~,sol] = ode45(@(t,x) dgl_d(t,x(1),x(2),c2,d2,m2,kopplung([1,5]),kopplung([2,6]),c3,d3,t0,h), [t_sol,t_sol + h], inital([3 4]));
    sys2 = sol(end,[1 2]);

end







%% Ideen
%{


Zeitpunkte der ode45 vorgeben - sol1 kann vorher initialisiert werden

Die anderen Übertragungsfunktionen implementieremn


Bei Force-Force müssen nur zwei u gespeichert werden, als initial
condition impelemtieren?

%}