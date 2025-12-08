clc
clear

load("OutputData\Log_1206_130320.mat")

h = Aero.Animation;

h.FramesPerSecond = 30;
h.TimeScaling = 1; 

idx1 = createBody(h, "C:\Users\emmat\Desktop\Year5\AERO626\Project\PlottingFunctions\Rocket Assembly.stl","");