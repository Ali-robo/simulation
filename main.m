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

n = 10000; %Anzahl Zeitschritte am Ende
h = 0.0001; %Zeitschritte

data = ff(n,h,sysPar,initial_conditions);

dataNumeric = calcNumericSol(sysPar,initial_conditions,linspace(0,h*n,n+1));

figure;

uNumeric = (sysPar.c3 .* (dataNumeric.x2-dataNumeric.x1) + sysPar.d3 .* (dataNumeric.v2 - dataNumeric.v1));
%plot(dataNumeric.time(2:n+1), data.uPre - data.uDebug(3:n+2)); hold on;
%plot([-h dataNumeric.time'],data.uDebug - uNumeric);hold off;
plot(dataNumeric.time, data.uDebug(2:n+2) - uNumeric); hold on;
plot(dataNumeric.time(2:n+1), data.uPre - uNumeric(2:n+1));


for va = ["x1","v1","x2","v2"]
    
    plot(dataNumeric.time, data.(va) - dataNumeric.(va)); hold on;
    legend(va);

end
legend("u","x1","v1","x2","v2");

% function lin = linInter(t,u0,u1,t0,h)
%     lin = u0 + (t-t0) * (u1-u0)/(h);
% end

%% plotting


figure;
time = linspace(0,h*n,n+1);

nexttile

title("System 1");

plot(time,data.x1 - dataNumeric.x1, "Color","r"); hold on;
plot(time,data.v1 - dataNumeric.v1, "Color","b");

%plot(dataNumeric.time, dataNumeric.x1, "Color","r","LineStyle","--");
%plot(dataNumeric.time, dataNumeric.v1, "Color","b","LineStyle","--");
hold off;

%ylim([-100,100]);
xlabel("time");
ylabel("x1, v1");
legend("x1", "v1", "x1x", "v1n");


nexttile

title("System 2");

plot(time,data.x2 - dataNumeric.x2, "Color","r"); hold on;
plot(time,data.v2 -dataNumeric.v2, "Color","b"); 

%plot(dataNumeric.time, dataNumeric.x2, "Color","r","LineStyle","--");
%plot(dataNumeric.time, dataNumeric.v2, "Color","b","LineStyle","--");
hold off;

%ylim([-100,100]);
xlabel("time");
ylabel("x2, v2");
legend("x2", "v2","x1x", "v1n");