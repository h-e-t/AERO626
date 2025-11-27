fileName = "C:\Users\emmat\Desktop\Year5\AERO626\Project\ThrustCurves\AeroTech_E24C.csv";

data = readmatrix(fileName);

tTime = data(:,1);
tThrust = data(:,2);

% Optionally ensure it starts at 0
if tTime(1) ~= 0
    tTime = [0; tTime];
    tThrust = [0; tThrust];
end

[path, name, ext] = fileparts(fileName); 
saveName = name + '.mat';

save(saveName, "tThrust", "tTime", '-mat'); 

