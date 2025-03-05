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
h_start = 1e-4; 
step = h_start;
count = 1000;

%dataNumeric = calcNumericSol(sysPar,initial_conditions,linspace(0,h_start*count*n,count*n+1));

dataFehler = zeros(count,2);
timeTotal = 0;

for i = 1:count

    tic

    dataCosim = ff(n,h_start*i,sysPar,initial_conditions);
    
    dataFehler(i,:) = localFehler(dataCosim,h_start*i,n,sysPar);
    
    timeTotal = timeTotal + toc;
    
    disp("Sim " + i + ", Time left " + (timeTotal/(i*60)*(count-i)));

end

 %% plot

figure;
ax = axes;

loglog(h_start*(1:count),dataFehler);
legend(ax,"x","v");
xlabel(ax,"h");
ylabel("error");
grid on;


%% functions

function fehler = globalFehler(dataCosim, dataNumeric,i)

    x1Numeric = dataNumeric.x1([linspace(1,i*n+1,((i*n)/i)+1)]);
    x2Numeric = dataNumeric.x2([linspace(1,i*n+1,((i*n)/i)+1)]);
    v1Numeric = dataNumeric.v1([linspace(1,i*n+1,((i*n)/i)+1)]);
    v2Numeric = dataNumeric.v2([linspace(1,i*n+1,((i*n)/i)+1)]);


    meanX = [sum(x1Numeric),sum(x2Numeric)]./n;
    meanV = [sum(v1Numeric),sum(v2Numeric)]./n;
    
    fehler(1) = calcFehler([x1Numeric,x2Numeric],[dataCosim.x1,dataCosim.x2],meanX);

    fehler(2) = calcFehler([v1Numeric,v2Numeric],[dataCosim.v1,dataCosim.v2],meanV);

end



function plotting(data,dataNumeric,h_start,n)

    figure;
    time = linspace(0,h_start*n,n+1);
    
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

function fehler = localFehler(dataCosim, h, n,sysPar)

    c1 = sysPar.c1; c2 = sysPar.c2; c3 = sysPar.c3;
    d1 = sysPar.d1; d2 = sysPar.d2; d3 = sysPar.d3;
    m1 = sysPar.m1; m2 = sysPar.m2;

    ode = @(t,x) [x(2); (-c1 * x(1) -d1 * x(2) + c3*(x(3)-x(1)) + d3 * (x(4) - x(2)))/m1;
              x(4); (-c2 * x(3) -d2 * x(4) - c3*(x(3)-x(1)) - d3 * (x(4) - x(2)))/m2];

    options =  odeset(RelTol=1e-10,AbsTol=1e-12);
    
    dataLokal = zeros(n+1,4);
    dataLokal(1,:) = [dataCosim.x1(1),dataCosim.v1(1),dataCosim.x2(1),dataCosim.v2(1)];

    for T = 1:n

        [~,temp] = ode45(ode,[0 h], [dataCosim.x1(T),dataCosim.v1(T),dataCosim.x2(T),dataCosim.v2(T)], options);
        dataLokal(T+1,:) = temp(end,:);

    end

    xMean = sum(dataLokal([1 3]),1)./n;
    vMean = sum(dataLokal([2 4]),1)./n;



    fehler(1) = calcFehler(dataLokal(:,[1 3]),[dataCosim.x1,dataCosim.x2],xMean);

    fehler(2) = calcFehler(dataLokal(:,[2 4]),[dataCosim.v1,dataCosim.v2],vMean);

end

function fehler = calcFehler(xGenau, x,xMean)

    fehler = sqrt( ...
        sum((xGenau(:,1) - x(:,1)).^2) / sum((xGenau(:,1) - xMean(1)).^2) + ...
        sum((xGenau(:,2) - x(:,2)).^2) / sum((xGenau(:,2) - xMean(2)).^2));

end


