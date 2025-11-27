clear
clc
rng(626)

function xk = nonLinearStatePropNoNoise(xkm)
    xk = xkm - 0.01 * sin(xkm); 
end

function xk = nonLinearStateProp(xkm, Pww)
    xk = xkm - 0.01 * sin(xkm) + randn()*chol(Pww)'; 
end

function zk = nonLinearState2Measurement(xk)
    zk = .5 * sin(2*xk); 
end

function F_X = linearizedDynamics(x)
    F_X = 1 - .01*cos(x); 
end 

function H_X = linearizedMeasurements(x)
    H_X = cos(2*x); 
end

Pww = .01^2; 

Pvv = .02; 

randomTruthSignals = nan(1,500); 
measurementHistory = nan(1,500); 



%%
intialState = 1.5;
randomTruthSignals(1) = intialState; 
measurementHistory(1) = nonLinearState2Measurement(intialState);

for idx = 2:500
        randomTruthSignals(idx) = nonLinearStateProp(randomTruthSignals(idx-1), .01^2); 
        measurementHistory(idx) = nonLinearState2Measurement(randomTruthSignals(idx));
end

%%

figure(Name="True Trajectory vs Measurements")
plot(randomTruthSignals, DisplayName="True Trajectory", Color="black", LineWidth=2)
hold on
scatter(1:500, measurementHistory,Marker="x", MarkerFaceColor="flat",MarkerEdgeColor="red", DisplayName="Measurements")
legend(Interpreter= 'LaTeX',FontSize=11)

ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)

% figure2pdf("True Trajectory vs Measurements.pdf")

%% Part 2 EKF

means = 1.5;
means = reshape(means, 1,1,[]);

covariances = .15^2; 
covariances = reshape(covariances,1,1,[]);

weights = 1;


EKF_Implemented = ExtendedGaussianMixtureFilter(means, covariances, weights);

meanHistory = nan(2,500); 

covHistory  = nan(2,500); 

numItems = 500;
for idx = 1:numItems
    [meanVal, cov] = EKF_Implemented.getCondProps(); 

    meanHistory(1, idx) = meanVal; 
    covHistory(1, idx)  = cov; 

    Zk = measurementHistory(idx);

    EKF_Implemented.predictPDF_Extended(@nonLinearStatePropNoNoise, @linearizedDynamics, Pww);
    EKF_Implemented.updatePDF_Extended(@nonLinearState2Measurement, @linearizedMeasurements, Pvv,Zk); 

    [meanVal, cov] = EKF_Implemented.getCondProps(); 

    meanHistory(2, idx) = meanVal; 
    covHistory(2, idx)  = cov; 
end


meanHistory     = reshape(meanHistory, 1,1000); 
covHistory      = reshape(covHistory, 1,1000);
arbitraryTime   = reshape([1:500; 1:500], 1, 1000); 

figure(Name="Extended Kalman Filte")
tiledlayout(2,1);

nexttile
plot(arbitraryTime, reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="True State", Color="black", LineWidth=2)
hold on
plot(arbitraryTime,meanHistory, DisplayName="EKF Mean History", Color="#808080", LineWidth=2)
ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")

nexttile
plot(arbitraryTime, meanHistory-reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="EKF Mean Error",Color="#808080", LineWidth=2, LineStyle="-")
hold on
plot(arbitraryTime,sqrt(covHistory)*3, DisplayName="$\pm3\sigma$", Color="red", LineWidth=2)
plot(arbitraryTime,- sqrt(covHistory)*3,HandleVisibility="off", Color="red", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Error",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")
% 
figure2pdf("EKF vs True State with updates.pdf")

%% Part 3
UKF_Implemented = ExtendedGaussianMixtureFilter(means, covariances, weights, .5, 2, 0);

meanHistory = nan(2,500); 

covHistory  = nan(2,500); 

numItems = 500;
for idx = 1:numItems
    [meanVal, cov] = UKF_Implemented.getCondProps(); 

    meanHistory(1, idx) = meanVal; 
    covHistory(1, idx)  = cov; 

    Zk = measurementHistory(idx);

    UKF_Implemented.predictPDF_Unscented(@nonLinearStatePropNoNoise, Pww);
    UKF_Implemented.updatePDF_Unscented(@nonLinearState2Measurement, Pvv,Zk); 

    [meanVal, cov] = UKF_Implemented.getCondProps(); 

    meanHistory(2, idx) = meanVal; 
    covHistory(2, idx)  = cov; 
end


meanHistory     = reshape(meanHistory, 1,1000); 
covHistory      = reshape(covHistory, 1,1000);
arbitraryTime   = reshape([1:500; 1:500], 1, 1000); 

figure(Name="Unscented Kalman Filter")
tiledlayout(2,1);

nexttile
plot(arbitraryTime, reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="True State", Color="black", LineWidth=2)
hold on
plot(arbitraryTime,meanHistory, DisplayName="UKF Mean History", Color="#808080", LineWidth=2)
ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")

nexttile
plot(arbitraryTime, meanHistory-reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="UKF Mean Error",Color="#808080", LineWidth=2, LineStyle="-")
hold on
plot(arbitraryTime,sqrt(covHistory)*3, DisplayName="$\pm3\sigma$", Color="red", LineWidth=2)
plot(arbitraryTime,- sqrt(covHistory)*3,HandleVisibility="off", Color="red", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Error",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")
% 
figure2pdf("UKF vs True State with updates.pdf")