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

%n = 1000; %Anzahl Zeitschritte am Ende
%h = 0.001; %Zeitschritte

plottingLocalFehler(dataFehler,0.00005*(1:4332))

%plottingLocalFehler(dataLocalFehler,linspace(h,h + (2 * 0.00005),3))

%ergebnisMain = fehlerU_h(sysPar,initial_conditions,h,n,1000);
 %data = fd(n,h,sysPar,initial_conditions,1);
 %dataNumeric = calcNumericSol(sysPar,initial_conditions,linspace(0,h*n,n+1));
% 
% 
%plotting(data,h,n);
%
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

function dataFehler = localFehler(sysPar, n, h_start,h_step, init)

    reps = 2;

    H = linspace(h_start,h_start + (reps * h_step),reps+1);


    dataFehler = struct();

    index =  1;
    time = 0;
    tic
    for h = H
        for k = 0:3
            
            dataFehler.ff.("k" + k)(index,:) = localFehler(ff(n,h,sysPar,init,k),h,n,sysPar);
            % dataFehler.df.("k" + k)(index,:) = localFehler(df(n,h,sysPar,init,k),h,n,sysPar);
            % dataFehler.fd.("k" + k)(index,:) = localFehler(fd(n,h,sysPar,init,k),h,n,sysPar);
            % dataFehler.dd.("k" + k)(index,:) = localFehler(dd(n,h,sysPar,init,k),h,n,sysPar);

        end

        time = time + toc;
        disp("noch " + (time)/((index)*60)*(reps+1-index) + " min");
       
        index = index +1;
        dataFehler.H(end) = h;
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
    
        fehler = zeros(1,2);
    
        fehler(1) = calcFehler(dataLokal(:,[1 3]),[dataCosim.x1,dataCosim.x2],xMean);
    
        fehler(2) = calcFehler(dataLokal(:,[2 4]),[dataCosim.v1,dataCosim.v2],vMean);
    
    end


    function fehler = calcFehler(xGenau, x,xMean)

    fehler = sqrt( ...
        sum((xGenau(:,1) - x(:,1)).^2) / sum((xGenau(:,1) - xMean(1)).^2) + ...
        sum((xGenau(:,2) - x(:,2)).^2) / sum((xGenau(:,2) - xMean(2)).^2));

    end



end

function plottingLocalFehler(dataFehler,H)

    index = 1;
    for koppl = fieldnames(dataFehler)'
        figure
        ax1 = subplot(1,2,1);
        ax2 = subplot(1,2,2);
        for k = fieldnames(dataFehler.(koppl{1}))'

            loglog(ax1,H,dataFehler.(koppl{1}).(k{1})(:,1));
            loglog(ax2,H,dataFehler.(koppl{1}).(k{1})(:,2));

            hold(ax1,"on");
            hold(ax2,"on");

        end
        
        % xlim(ax1,[1e-5 , 1e-2])
        % xlim(ax2,[1e-5 , 1e-2])

        hold(ax1,"off");
        hold(ax2,"off");

        grid(ax1,"minor");
        grid(ax2,"minor");
        
        title(ax1,"locError x, Kopplung: " + (koppl{1}));
        title(ax2,"locError v, Kopplung: " + (koppl{1}));

        legend(ax1,fieldnames(dataFehler.(koppl{1})))
        legend(ax2,fieldnames(dataFehler.(koppl{1})))

        index = index + 2;
    end

end