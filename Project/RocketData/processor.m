cnAlpha = readmatrix("CNalpha.csv"); 
copData = readmatrix("COPvsAOA.csv"); 

aoa = [-cnAlpha(end:-1:2,1);cnAlpha(:,1)]; 
cna = [cnAlpha(end:-1:2,2);cnAlpha(:,2)]; % conversion to kg 
cop = [copData(end:-1:2,2);copData(:,2)]/100; % conversion to meters

save("vehicleData.mat", "aoa", "cop", "cna"); 

