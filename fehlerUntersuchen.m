%fehler Untersuchen


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

n = 10000; 
h_start = 1e-4; %Zeitschritte
step = h_start;
count = 10;


dataFehler = struct("h", zeros(count,1),"x",zeros(count,1),"v",zeros(count,1),"test",zeros(count,1));

dataNumeric = calcNumericSol(sysPar,initial_conditions,linspace(0,h*count*n,count*n+1));

% meanX = [sum([dataNumeric.x1,dataNumeric.x2],1)]./n;
% meanV = [sum([dataNumeric.v1,dataNumeric.v2],1)]./n;

for i = 1:count

    dataCosim = ff_const(n,h*i,sysPar,initial_conditions);

    x1Numeric = dataNumeric.x1([linspace(1,i*n+1,((i*n)/i)+1)]);
    x2Numeric = dataNumeric.x2([linspace(1,i*n+1,((i*n)/i)+1)]);
    v1Numeric = dataNumeric.v1([linspace(1,i*n+1,((i*n)/i)+1)]);
    v2Numeric = dataNumeric.v2([linspace(1,i*n+1,((i*n)/i)+1)]);

    plotting(dataCosim, struct("x1", x1Numeric,"v1", v1Numeric,"x2", x2Numeric, "v2", v2Numeric,"time", linspace(0,h*n*i,n+1)),h*i,n);

    meanX = [sum(x1Numeric),sum(x2Numeric)]./n;
    meanV = [sum(v1Numeric),sum(v2Numeric)]./n;
    
    dataFehler.x(i) = globalFehler([x1Numeric,x2Numeric],[dataCosim.x1,dataCosim.x2],meanX);

    dataFehler.v(i) = globalFehler([v1Numeric,v2Numeric],[dataCosim.v1,dataCosim.v2],meanV);

    dataFehler.test(i) = abs(mean(dataCosim.x1 - x1Numeric));
    
    disp("Sim " + i);
    
end



%% functions


function fehler = globalFehler(xGenau, x,xMean)

    fehler = sqrt( ...
        sum((xGenau(:,1) - x(:,1)).^2) / sum((xGenau(:,1) - xMean(1)).^2) + ...
        sum((xGenau(:,2) - x(:,2)).^2) / sum((xGenau(:,2) - xMean(2)).^2));

end

function plotting(data,dataNumeric,h,n)

    figure;
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

end





 %sum(sum((xGenau - x).^2,1) / sum((xGenau - xMean).^2,1));

 %% plot

figure;
ax = axes;

loglog(h*(1:count),dataFehler.x); hold on;
loglog(h*(1:count),dataFehler.v); 
loglog(h*(1:count),dataFehler.test); hold off;
legend(ax,"x","v","absoluter Fehler x1 mean");
xlabel(ax,"h");
ylabel("global error");
grid on;
