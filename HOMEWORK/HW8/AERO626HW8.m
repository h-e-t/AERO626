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




initialState = 1.5;
randomTruthSignals(1) = initialState; 
measurementHistory(1) = nonLinearState2Measurement(initialState);

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

tic
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
EKFTTC = toc;

EKFError =  meanHistory(2,:)-randomTruthSignals(1,:); 

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
tic
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
UKFTTC = toc;
UKFError =  meanHistory(2,:)-randomTruthSignals(1,:); 

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

%% Part 4
initDist = ExtendedGaussianMixtureFilter(1.5, .15^2, 1);
particleFilter = BPFilter(10000, initDist);

meanHistory = nan(2,500); 

covHistory  = nan(2,500); 

NeffHistory = nan(2,500); 

tic
numItems = 500;
for idx = 1:numItems
    [meanVal, cov] = particleFilter.getProps(); 
    
    meanHistory(1, idx) = meanVal; 
    covHistory(1, idx)  = cov; 
    NeffHistory(1,idx)  = particleFilter.calculateEffectiveParticles();

    Zk = measurementHistory(idx);

    particleFilter.propagate(@nonLinearStateProp, Pww);
    particleFilter.update(@nonLinearState2Measurement, Pvv,Zk); 

    [meanVal, cov] = particleFilter.getProps(); 

    meanHistory(2, idx) = meanVal; 
    covHistory(2, idx)  = cov; 
    NeffHistory(2,idx)  = particleFilter.calculateEffectiveParticles();
end
BPF_NR_TTC = toc;

BPFNotResampledError =  meanHistory(2,:)-randomTruthSignals(1,:); 

meanHistory     = reshape(meanHistory, 1,1000); 
covHistory      = reshape(covHistory, 1,1000);
NeffHistory     = reshape(NeffHistory, 1,1000);

arbitraryTime   = reshape([1:500; 1:500], 1, 1000); 

figure(Name="Bootstrap Particle Filter")
tiledlayout(2,1);

nexttile
plot(arbitraryTime, reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="True State", Color="black", LineWidth=2)
hold on
plot(arbitraryTime,meanHistory, DisplayName="BPF Mean History", Color="#808080", LineWidth=2)
ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")
nexttile
plot(arbitraryTime, meanHistory-reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="BPF Mean Error",Color="#808080", LineWidth=2, LineStyle="-")
hold on
plot(arbitraryTime,sqrt(covHistory)*3, DisplayName="$\pm3\sigma$", Color="red", LineWidth=2)
plot(arbitraryTime,- sqrt(covHistory)*3,HandleVisibility="off", Color="red", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Error",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")

figure2pdf("BPF vs True State No Resampling.pdf")


figure(Name="Neff over time")

plot(arbitraryTime,NeffHistory, DisplayName="$N_{eff}$", Color="#808080", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Number of Particles",Interpreter= 'LaTeX',FontSize=11)

legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")


figure2pdf("Neff over Time No Resampling.pdf")


%% Part 5
initDist = ExtendedGaussianMixtureFilter(1.5, .15^2, 1);
particleFilter = BPFilter(10000, initDist, true, 1000);

meanHistory = nan(2,500); 

covHistory  = nan(2,500); 

NeffHistory = nan(2,500); 
tic
numItems = 500;
for idx = 1:numItems
    [meanVal, cov] = particleFilter.getProps(); 
    
    meanHistory(1, idx) = meanVal; 
    covHistory(1, idx)  = cov; 
    NeffHistory(1,idx)  = particleFilter.calculateEffectiveParticles();

    Zk = measurementHistory(idx);

    particleFilter.propagate(@nonLinearStateProp, Pww);
    particleFilter.update(@nonLinearState2Measurement, Pvv,Zk); 

    [meanVal, cov] = particleFilter.getProps(); 

    meanHistory(2, idx) = meanVal; 
    covHistory(2, idx)  = cov; 
    NeffHistory(2,idx)  = particleFilter.calculateEffectiveParticles();
end
BPF_R_TTC = toc;

numTimesResampled = particleFilter.resampleCount

BPFResampledError =  meanHistory(2,:)-randomTruthSignals(1,:); 

meanHistory     = reshape(meanHistory, 1,1000); 
covHistory      = reshape(covHistory, 1,1000);
NeffHistory     = reshape(NeffHistory, 1,1000);

arbitraryTime   = reshape([1:500; 1:500], 1, 1000); 

figure(Name="Bootstrap Particle Filter")
tiledlayout(2,1);

nexttile
plot(arbitraryTime, reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="True State", Color="black", LineWidth=2)
hold on
plot(arbitraryTime,meanHistory, DisplayName="BPF Mean History", Color="#808080", LineWidth=2)
ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")
nexttile
plot(arbitraryTime,meanHistory-reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="BPF Mean Error",Color="#808080", LineWidth=2, LineStyle="-")
hold on
plot(arbitraryTime,sqrt(covHistory)*3, DisplayName="$\pm3\sigma$", Color="red", LineWidth=2)
plot(arbitraryTime,- sqrt(covHistory)*3,HandleVisibility="off", Color="red", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Error",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")
figure2pdf("BPF vs True State With Resampling.pdf")

figure(Name="Neff over time")

plot(arbitraryTime,NeffHistory, DisplayName="$N_{eff}$", Color="#808080", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Number of Particles",Interpreter= 'LaTeX',FontSize=11)

legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")

figure2pdf("Neff over Time With Resampling.pdf")

%% Part 6
ekfrmse = rms(EKFError);
ukfrmse = rms(UKFError);
bpfNRrmse = rms(BPFNotResampledError);
bpfRrmse = rms(BPFResampledError);

RMSErrors = [ekfrmse; ukfrmse; bpfNRrmse; bpfRrmse];

ekfmae = mean(abs(EKFError));
ukfmae = mean(abs(UKFError));
bpfNRmae = mean(abs(BPFNotResampledError));
bpfRmae = mean(abs(BPFResampledError));

MAErrors = [ekfmae; ukfmae; bpfNRmae; bpfRmae];

Algorithm=["EKF" "UKF" "BPF" "Resampled BPF"]';
TimesToCompute = [EKFTTC, UKFTTC, BPF_NR_TTC, BPF_R_TTC]';
T = table(Algorithm, RMSErrors, MAErrors, TimesToCompute)

%% Part 7

EKF_MC_RMSE = nan(1,500); 
EKF_MC_MAE = nan(1,500); 

UKF_MC_RMSE = nan(1,500); 
UKF_MC_MAE = nan(1,500); 

BPF_MC_RMSE = nan(1,500); 
BPF_MC_MAE = nan(1,500); 

initDist = ExtendedGaussianMixtureFilter(1.5, .15^2, 1);

d1 = 500; 
d2 = 500; 
randomTruthSignals = nan(d1,d2); 

for idx = 1:d1
    intialState = initDist.sampleDistribution();
    randomTruthSignals(idx, 1) = intialState; 
    for jdx=2:d2
        randomTruthSignals(idx, jdx) = nonLinearStateProp(randomTruthSignals(idx, jdx-1), .01^2); 
    end
end

%%

for run = 1:500
    EKF_Implemented = ExtendedGaussianMixtureFilter(means, covariances, weights);
    meanHistory = nan(1,500); 
    covHistory  = nan(1,500); 
    
    numItems = 500;
    for idx = 1:numItems
        [meanVal, cov] = EKF_Implemented.getCondProps(); 
    
        meanHistory(1, idx) = meanVal; 
        covHistory(1, idx)  = cov; 
    
        Zk =  nonLinearState2Measurement(randomTruthSignals(run, idx));
    
        EKF_Implemented.predictPDF_Extended(@nonLinearStatePropNoNoise, @linearizedDynamics, Pww);
        EKF_Implemented.updatePDF_Extended(@nonLinearState2Measurement, @linearizedMeasurements, Pvv,Zk); 
    end

    error = meanHistory-randomTruthSignals(1,:); 
    
    EKF_MC_RMSE(run) = rms(error); 
    EKF_MC_MAE(run)  = mean(abs(error)); 

end
%%

for run = 1:500
    UKF_Implemented = ExtendedGaussianMixtureFilter(means, covariances, weights, .5, 2, 0);
    
    meanHistory = nan(1,500); 
    
    covHistory  = nan(1,500); 
    
    numItems = 500;
    for idx = 1:numItems
        [meanVal, cov] = UKF_Implemented.getCondProps(); 
    
        meanHistory(1, idx) = meanVal; 
        covHistory(1, idx)  = cov; 
    
        Zk =  nonLinearState2Measurement(randomTruthSignals(run, idx));
    
        UKF_Implemented.predictPDF_Unscented(@nonLinearStatePropNoNoise, Pww);
        UKF_Implemented.updatePDF_Unscented(@nonLinearState2Measurement, Pvv,Zk); 
    end

    error = meanHistory-randomTruthSignals(1,:); 
    
    UKF_MC_RMSE(run) = rms(error); 
    UKF_MC_MAE(run)  = mean(abs(error)); 

end

%%
% timeTaken = 5.5; 
% 
% for run = 1:500
%     if mod(run,10) == 0
%         clc
%         disp("Trial " + run)
%         disp("Time remamining = " + ((500 - run) * timeTaken)/60 + " minutes ")
%     end
% 
%     initDist = ExtendedGaussianMixtureFilter(1.5, .15^2, 1);
%     particleFilter = BPFilter(10000, initDist, true, 1000);
% 
%     meanHistory = nan(1,500); 
% 
%     covHistory  = nan(1,500); 
%     tic
%     numItems = 500;
%     for idx = 1:numItems
%         [meanVal, cov] = particleFilter.getProps(); 
% 
%         meanHistory(1, idx) = meanVal; 
%         covHistory(1, idx)  = cov; 
% 
%         Zk = nonLinearState2Measurement(randomTruthSignals(run, idx));
% 
%         particleFilter.propagate(@nonLinearStateProp, Pww);
%         particleFilter.update(@nonLinearState2Measurement, Pvv,Zk); 
% 
%         [meanVal, cov] = particleFilter.getProps(); 
% 
%     end
%     timeTaken = toc;
% 
%     error = meanHistory-randomTruthSignals(1,:); 
% 
%     BPF_MC_RMSE(run) = rms(error); 
%     BPF_MC_MAE(run)  = mean(abs(error)); 
% end

%%
figure()

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

set(groot, 'defaultTextFontSize', 11);
set(groot, 'defaultAxesFontSize', 11)

data = [
    mean(EKF_MC_RMSE)  mean(EKF_MC_MAE);
    mean(UKF_MC_RMSE)  mean(UKF_MC_MAE);
    mean(BPF_MC_RMSE)  mean(BPF_MC_MAE)
];

figure;
b = bar(data, 0.8, 'grouped');
ylabel("Error")
legend("Root Mean Squared Error", "Mean Absolute Error", 'Interpreter', 'LaTeX')

b(1).FaceColor = "#7C0A02"; % RMSE
b(2).FaceColor = "#808080"; % MAE

grid on
xticklabels(["EKF" "UKF" "BPF"]);


figure2pdf("Filter MC Comparison.pdf")

%%

figure()

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

set(groot, 'defaultTextFontSize', 11);
set(groot, 'defaultAxesFontSize', 11)

subplot(3,2,1)
histogram(EKF_MC_RMSE, 100,FaceColor = "#7C0A02")
title("EKF RMSE")

subplot(3,2,2)
histogram(EKF_MC_MAE, 100,FaceColor = "#7C0A02")
title("EKF MAE")

subplot(3,2,3)
histogram(UKF_MC_RMSE, 100,FaceColor = "#7C0A02")
title("UKF RMSE")

subplot(3,2,4)
histogram(UKF_MC_MAE, 100,FaceColor = "#7C0A02")
title("UKF MAE")

subplot(3,2,5)
histogram(BPF_MC_RMSE, 100,FaceColor = "#7C0A02")
title("BPF RMSE")

subplot(3,2,6)
histogram(BPF_MC_MAE, 100,FaceColor = "#7C0A02")
title("BPF MAE")

figure2pdf("Complete Filter MC Comparison.pdf")
