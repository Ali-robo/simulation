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
h = 0.001; %Zeitschritte


%ergebnisMain = fehlerU_h(sysPar,initial_conditions,h,n,1000);


 data = fd(n,h,sysPar,initial_conditions,1);
 dataNumeric = calcNumericSol(sysPar,initial_conditions,linspace(0,h*n,n+1));
% 
% 
plotting(data,dataNumeric,h,n);

%fehlerU_h(sysPar,data,dataNumeric,n,h);
%fehlerU(sysPar,data,dataNumeric,n,h);

%% plotting

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


function fehlerU(sysPar,data,dataNumeric,n,h)

    figure;

    ax = axes;
    uNumeric = (sysPar.c3 .* (dataNumeric.x2-dataNumeric.x1) + sysPar.d3 .* (dataNumeric.v2 - dataNumeric.v1));
    plot(ax,dataNumeric.time, data.uDebug(2:n+2) - uNumeric); hold on;
    plot(ax,dataNumeric.time(2:n+1), data.uPre - uNumeric(2:n+1));
    
    
    for va = ["x1","v1","x2","v2"]
        
        plot(ax,dataNumeric.time, data.(va) - dataNumeric.(va)); hold on;
      
    end

    legend(ax,"u","uPre","x1","v1","x2","v2");
    title(ax,"Absoluter Fehler von den Systemen und der Kopplung über die Zeit","h: " + h);

end


function value = fehlerU_h(sysPar,init,h_start,n,count)

    fehlerUmax = zeros(count,1);
    fehlerUmean = zeros(count,1);

    dataNumeric = calcNumericSol(sysPar,init,linspace(0,h_start*count*n,count*n+1));

    uTemp = (sysPar.c3 .* (dataNumeric.x2-dataNumeric.x1) + sysPar.d3 .* (dataNumeric.v2 - dataNumeric.v1));

    uNumeric = zeros(n+1,count);

    uSim = zeros(n+1,count);

    for i = 1:count

        uNumeric(:,i) = uTemp([linspace(1,i*n+1,((i*n)/i)+1)]);

    end

    parfor i = 1:count

        uSim(:,i) = ff(n,h_start*i,sysPar,init).uDebug(2:n+2);

        disp("Sim" + i);

    end


    for i = 1:count

        fehlerUmax(i) = max(abs(uSim(:,i) - uNumeric(:,i)));
        fehlerUmean(i) = mean(abs(uSim(:,i) - uNumeric(:,i)));

    end

    plot(h_start*(1:count),fehlerUmax); hold on;
    plot(h_start*(1:count), fehlerUmean); hold off;

    value = struct("max" , fehlerUmax, "mean", fehlerUmean);


end