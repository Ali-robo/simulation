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

n = 1000; %Anzahl Zeitschritte am Ende
h = 0.0001; %Zeitschritte



fehlerU_h(h,n,sysPar,initial_conditions,100);

%data = ff(n,h,sysPar,initial_conditions);

%dataNumeric = calcNumericSol(sysPar,initial_conditions,linspace(0,h*n,n+1));


legend("u","x1","v1","x2","v2");


%% plotting

function plotting(data,h,n)

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


function fehlerU(data,dataNumeric,n)

    figure;

    uNumeric = (sysPar.c3 .* (dataNumeric.x2-dataNumeric.x1) + sysPar.d3 .* (dataNumeric.v2 - dataNumeric.v1));
    plot(dataNumeric.time, data.uDebug(2:n+2) - uNumeric); hold on;
    plot(dataNumeric.time(2:n+1), data.uPre - uNumeric(2:n+1));
    
    
    for va = ["x1","v1","x2","v2"]
        
        plot(dataNumeric.time, data.(va) - dataNumeric.(va)); hold on;
      
    end

     legend("u","uPre","x1","v1","x2","v2");

end



function fehlerU_h(h_start,n,sysPar,init,count)

    fehlerUmax = zeros(count,1);
    fehlerUmean = zeros(count,1);

    dataNumeric = calcNumericSol(sysPar,init,linspace(0,h_start*count*n,count*n+1));

    uTemp = (sysPar.c3 .* (dataNumeric.x2-dataNumeric.x1) + sysPar.d3 .* (dataNumeric.v2 - dataNumeric.v1));


    for i = 1:count

        uNumeric(:,i) = uTemp([linspace(1,i*n+1,((i*n)/i)+1)]);

    end

    parfor i = 1:count

        data = ff(n,h_start*i,sysPar,init);

        fehlerUmax(i) = max(abs(data.uDebug(2:n+2) - uNumeric(:,i)));
        fehlerUmean(i) = mean(abs(data.uDebug(2:n+2) - uNumeric(:,i)));

        disp("Sim" + i);

    end

    plot(h_start*(1:count),fehlerUmax); hold on;
    plot(h_start*(1:count), fehlerUmean); hold off;


end