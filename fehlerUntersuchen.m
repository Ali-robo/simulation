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

n = 1000; 
h = 1e-5; %Zeitschritte
count = 15000;


dataFehler = struct("mean",struct);

dataNumeric = calcNumericSol(sysPar,initial_conditions,linspace(0,h*count*n,count*n+1));

for i = 1:count

    dataCosim = df(n,h*i,sysPar,initial_conditions);

    for va = ["x1","v1","x2","v2"]

         dataFehler.("Sim" +i).("h") = h*i;

         dataFehler.("Sim" + i).(va) = abs(dataCosim.(va) - dataNumeric.(va)([linspace(1,i*n+1,((i*n)/i)+1)]));

         dataFehler.mean.(va)(i) = mean(dataFehler.("Sim" + i).(va)); % mean ergibt glaube ich kein Sinn
         dataFehler.max.(va)(i) = max(dataFehler.("Sim" + i).(va));

    end
    
    disp("Sim " + i);
    
end

%% plot

figure;
ax(1) = subplot(1,2,1);
ax(2) = subplot(1,2,2);

subplot(ax(1))
loglog(h*(1:count),dataFehler.mean.x1,"--r"); hold on;
loglog(h*(1:count),dataFehler.mean.x2,"--b")
loglog(h*(1:count),dataFehler.mean.v1,"-r")
loglog(h*(1:count),dataFehler.mean.v2,"-b"); hold off;
legend("x1","v1","x2","v2");
ax(1).Title.String = "Durschnittlicher Fehler über Zeitschritte";
xlabel("h");
ylabel("Durschnittlicher Fehler");
grid on;


subplot(ax(2));
loglog(h*(1:count),dataFehler.max.x1,"--r"); hold on;
loglog(h*(1:count),dataFehler.max.x2,"--b")
loglog(h*(1:count),dataFehler.max.v1,"-r")
loglog(h*(1:count),dataFehler.max.v2,"-b"); hold off;
legend("x1","v1","x2","v2");
title(ax(2),"Maximaler Fehler über Zeitschritte");
xlabel("h");
ylabel("Maximaler Fehler");
grid on;