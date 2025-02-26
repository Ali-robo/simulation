sysPar = struct( ...
    'c1', 400, ...
    'c2', 300, ...
    'c3', 80, ...
    'd1', 0.2, ...
    'd2', 0.3, ...
    'd3', 15, ...
    'm1', 6, ...
    'm2',6 ...
    );

initial_conditions = [4, 2, 1, 3];

n = 10; %Anzahl Zeitschritte am Ende
h = 0.1; %Zeitschritte

data = ff(n,h,sysPar,initial_conditions);

dataNumeric = calcNumericSol(sysPar,initial_conditions,h,n);


%% plotting

time = linspace(0,h*n,n+1);

nexttile

title("System 1");

plot(time,data.x1, "Color","r"); hold on;
plot(time,data.v1, "Color","b");

plot(dataNumeric.time, dataNumeric.x1, "Color","r","LineStyle","--");
plot(dataNumeric.time, dataNumeric.v1, "Color","b","LineStyle","--");
hold off;

ylim([-100,100]);
xlabel("time");
ylabel("x1, v1");
legend("x1", "v1", "x1x", "v1n");


nexttile

title("System 2");

plot(time,data.x2, "Color","r"); hold on;
plot(time,data.v2, "Color","b"); 

plot(dataNumeric.time, dataNumeric.x2, "Color","r","LineStyle","--");
plot(dataNumeric.time, dataNumeric.v2, "Color","b","LineStyle","--");
hold off;

ylim([-100,100]);
xlabel("time");
ylabel("x2, v2");
legend("x2", "v2","x1x", "v1n");