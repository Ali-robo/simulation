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

n = 1000; %Anzahl Zeitschritte am Ende
h = 0.000001; %Zeitschritte
count = 1000;


dataFehler = struct("mean",struct);

for i = 1:count

    dataCosim = ff(n,h*i,sysPar,initial_conditions);
    
    dataNumeric = calcNumericSol(sysPar,initial_conditions,h*i,n);

    
    
    for va = ["x1","v1","x2","v2"]

         dataFehler.("Sim" +i).("h") = h*i;

         dataFehler.("Sim" + i).(va) = abs(dataCosim.(va) - dataNumeric.(va));

         dataFehler.mean.(va)(i) = mean(dataFehler.("Sim" + i).(va));

    end
    
    disp("Sim " + i);
    

end


loglog(h*(1:count),dataFehler.mean.x1,"--r"); hold on;
loglog(h*(1:count),dataFehler.mean.x2,"--b")
loglog(h*(1:count),dataFehler.mean.v1,"-r")
loglog(h*(1:count),dataFehler.mean.v2,"-b"); hold off;
grid on;