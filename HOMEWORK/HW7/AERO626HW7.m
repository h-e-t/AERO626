clear
clc

function nextStep = nonLinearStatePropNoNoise(xkm)
    nextStep = xkm - 0.01 * sin(xkm); 
end

function nextStep = nonLinearStateProp(xkm, Pww)
    nextStep = xkm - 0.01 * sin(xkm) + randn()*chol(Pww)'; 
end

function meanMeas = nonLinearState2Measurement(xk)
    meanMeas = .5 * sin(2*xk); 
end

function F_X = linearizedDynamics(x)
    F_X = 1 - .01*cos(x); 
end 

function H_X = linearizedMeasurements(x)
    H_X = cos(2*x); 
end

Pww = .01^2; 

Pvv = .02; 

means = [1.0909 1.1818 1.2727 1.3636 1.4545 1.5455 1.6364 1.7273 1.8182 1.9091];
means = reshape(means, 1,1,[]);

covariances = ones(1,10)*.0696^2; 
covariances = reshape(covariances,1,1,[]);

weights = ones(1,10)*.1;

initialDistribution=ExtendedGaussianMixtureFilter(means, covariances, weights);

% dist=initialDistribution.distribution(); 
% 
% plot(dist.xEval, dist.distributions)
%%


% Part 1

d1 = 1; 
d2 = 500; 
randomTruthSignals = nan(d1,d2); 

for idx = 1:d1
    intialState = initialDistribution.sampleDistribution();
    randomTruthSignals(idx, 1) = intialState; 
    for jdx=2:d2
        randomTruthSignals(idx, jdx) = nonLinearStateProp(randomTruthSignals(idx, jdx-1), .01^2); 
    end
end

%%



figure(name="MC Sim");
tiledlayout(2,1); 
ax1 = nexttile;
ax2 = nexttile; 

scatter(ax1, 1:d2, randomTruthSignals(1,:), .2 ,"red", DisplayName="Monte Carlo Truth", Marker="x")
hold(ax1, "on")

for idx = 2:d1
    scatter(ax1, 1:d2, randomTruthSignals(idx,:),.1, "red", "HandleVisibility","off")
end

meanState = mean(randomTruthSignals, 1);
stateVariance = std(randomTruthSignals, 1).^2;

plot(ax1, 1:d2, meanState, LineWidth=2, Color='k', DisplayName="Sample Mean")
plot(ax2, 1:d2, stateVariance, LineWidth=2, Color='k', HandleVisibility="off")

legend(ax1, "show",Interpreter= 'LaTeX',FontSize=11)
xlabel(ax2, "Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel(ax1, "Value",Interpreter= 'LaTeX',FontSize=11)
ylabel(ax2, "Variance",Interpreter= 'LaTeX',FontSize=11)

figure2pdf("MC_Sim.pdf")

%% Part 2

Part2Distribution = ExtendedGaussianMixtureFilter(means, covariances, weights);

meanHistory = nan(1,500); 
covHistory  = nan(1,500); 
measurementHistory = nan(1,500); 

distributionHistories = nan(500,length(Part2Distribution.distribution().xEval)); 
evalHistories = nan(500, length(Part2Distribution.distribution().xEval));

figure(Name="GMPDF Plot Over Time")
tiledlayout()
ax = gca; 
hold(ax, "on")

numItems = 500;
for idx = 1:numItems
    [meanVal, cov] = Part2Distribution.getCondProps(); 

    meanHistory(idx) = meanVal; 
    covHistory(idx)  = cov; 
    measurementHistory(idx) = nonLinearState2Measurement(randomTruthSignals(1,idx));

    Part2Distribution.predictPDF_Extended(@nonLinearStatePropNoNoise, @linearizedDynamics,Pww);


    if any(idx==[1 0:25:500])
        dist = Part2Distribution.distribution();  
        distributionHistories(idx, :) = dist.distributions(1,:);
        evalHistories(idx, :) = dist.xEval;

        if any(idx == [1 250 500])
            plot3(ax, dist.xEval(1,:), idx*ones(size(dist.xEval(1,:))), dist.distributions(1,:), "Color",'r',DisplayName="EGM PDF "+idx)
        else
            plot3(ax, dist.xEval(1,:), idx*ones(size(dist.xEval(1,:))), dist.distributions(1,:), "Color",'b',DisplayName="EGM PDF",HandleVisibility="off")
        end

    end
end
legend(ax, "show",Interpreter= 'LaTeX',FontSize=11)
grid(ax)
plot3(ax, meanHistory(1:numItems), 1:numItems, zeros(1,numItems), DisplayName="Monte Carlo Mean State", Color="black",LineWidth=1.5)
plot3(ax, meanState(1,1:numItems), 1:numItems, .01*ones(1,numItems), DisplayName="EGMF Mean State", Color="#888888",LineWidth=3, LineStyle="--")
xlabel("$x$",Interpreter= 'LaTeX',FontSize=11)
ylabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
zlabel("$p(x)$",Interpreter= 'LaTeX',FontSize=11)
view(ax, 25,60)

figure2pdf("GMPDF Over Time.pdf")
%%

figure(Name="Extended GMF No Updates")
tiledlayout(2,1);

nexttile
plot(meanState, DisplayName="True State", Color="black", LineWidth=2)
hold on
plot(meanHistory, DisplayName="UGMF Mean History", Color="#808080", LineWidth=2)
ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")

nexttile
plot(meanHistory-meanState, DisplayName="UGMF Mean Error",Color="#808080", LineWidth=2, LineStyle="-")
hold on
plot(sqrt(covHistory)*3, DisplayName="$\pm3\sigma$", Color="red", LineWidth=2)
plot(- sqrt(covHistory)*3,HandleVisibility="off", Color="red", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Error",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")

figure2pdf("EMPDF Over Time No Updates.pdf")



%% Part 3
figure(Name = "Mean n Variance over Time")
tiledlayout(2,1)

nexttile
plot(meanHistory, DisplayName="EGMF Mean", Color="black", LineWidth=2); 
hold on;
plot(meanState, DisplayName="Monte Carlo Mean", Color="red", LineStyle="--", LineWidth=2)
ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)

nexttile
plot(stateVariance, Color="black", LineWidth=2,DisplayName="Monte Carlo Variance")
hold on 
ylim([0,.080])
plot(covHistory, Color="red", LineStyle="--", LineWidth=2,DisplayName="Extended GMF Variance")

ylabel("Variance",Interpreter= 'LaTeX',FontSize=11)
xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)

updateOnlyEGMFrmse = rms(meanHistory-meanState);
figure2pdf("EGMF Mean and Variance vs MC and Variance no updates.pdf")

%%
figure(Name="True Trajectory vs Measurements")
plot(randomTruthSignals(1,:), DisplayName="True Trajectory", Color="black", LineWidth=2)
hold on
scatter(1:500, measurementHistory,Marker="x", MarkerFaceColor="flat",MarkerEdgeColor="red", DisplayName="Measurements")
legend(Interpreter= 'LaTeX',FontSize=11)

ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)

figure2pdf("True Trajectory vs Measurements.pdf")


%% Part 4

Part4Distribution = ExtendedGaussianMixtureFilter(means, covariances, weights);

meanHistory = nan(1,500); 
covHistory  = nan(1,500); 
measurementHistory = nan(1,500); 

distributionHistories = nan(500,length(Part4Distribution.distribution().xEval)); 
evalHistories = nan(500, length(Part4Distribution.distribution().xEval));

figure(Name="GMPDF Plot Over Time")
ax = gca; 

hold(ax, "on")

numItems = 500;
for idx = 1:numItems
    [meanVal, cov] = Part4Distribution.getCondProps(); 

    meanHistory(idx) = meanVal; 
    covHistory(idx)  = cov; 

    Zk = nonLinearState2Measurement(randomTruthSignals(1,idx));

    Part4Distribution.predictPDF_Extended(@nonLinearStatePropNoNoise, @linearizedDynamics,Pww);
    Part4Distribution.updatePDF_Extended(@nonLinearState2Measurement, @linearizedMeasurements, Pvv,Zk); 


    if any(idx==1:25:500)
        dist = Part4Distribution.distribution();  
        distributionHistories(idx, :) = dist.distributions(1,:);
        evalHistories(idx, :) = dist.xEval;

        if (idx == 1)
            plot3(ax, dist.xEval(1,:), idx*ones(size(dist.xEval(1,:))), dist.distributions(1,:), "Color",'b',DisplayName="GM PDF")
        else
            plot3(ax, dist.xEval(1,:), idx*ones(size(dist.xEval(1,:))), dist.distributions(1,:), "Color",'b',DisplayName="GM PDF",HandleVisibility="off")
        end

    end
end
grid(ax, "on")
plot3(ax, meanHistory, 1:numItems, zeros(1,numItems), DisplayName="EGMF Mean History", Color="#808080", LineWidth=1.5, LineStyle="-")

plot3(ax, meanHistory + sqrt(covHistory)*3, 1:numItems, zeros(1,numItems), DisplayName="$\pm3\sigma$", Color="red", LineWidth=1.5)
plot3(ax, meanHistory - sqrt(covHistory)*3, 1:numItems, zeros(1,numItems), HandleVisibility="off", Color="red", LineWidth=1.5)

plot3(ax, randomTruthSignals(1,:), 1:numItems, zeros(1,numItems), DisplayName="True State", Color="black", LineWidth=1.5)

xlabel("$x$",Interpreter= 'LaTeX',FontSize=11)
ylabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
zlabel("$p(x)$",Interpreter= 'LaTeX',FontSize=11)
legend(ax,Interpreter= 'LaTeX',FontSize=11)

view(ax, 25,60)

figure2pdf("EGMF vs True State pdf with updates.pdf")
%%
figure(Name ="EGMF vs true state")
tiledlayout(2,1);

nexttile
plot(randomTruthSignals(1,:), DisplayName="True State", Color="black", LineWidth=2)
hold on
plot(meanHistory, DisplayName="EGMF Mean History", Color="#808080", LineWidth=2)
legend(Interpreter= 'LaTeX',FontSize=11)
ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
grid("on")

nexttile
plot(meanHistory-randomTruthSignals(1,:), DisplayName="EGMF Mean Error",Color="#808080", LineWidth=2, LineStyle="-")
hold on
plot(sqrt(covHistory)*3, DisplayName="$\pm3\sigma$", Color="red", LineWidth=2)
plot(- sqrt(covHistory)*3,HandleVisibility="off", Color="red", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Error",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")

figure2pdf("EGMF vs True State mean and variance with updates.pdf")

EGMF_rmse = rms(meanHistory-randomTruthSignals(1,:))


%% Part 5

means = 1.5;
means = reshape(means, 1,1,[]);

covariances = .15^2; 
covariances = reshape(covariances,1,1,[]);

weights = 1;


Part5Distribution = ExtendedGaussianMixtureFilter(means, covariances, weights, .5, 2, 0);
Part5Distribution.predictPDF_Unscented(@nonLinearStatePropNoNoise, Pww)

meanHistory = nan(2,500); 

covHistory  = nan(2,500); 

measurementHistory = nan(1,500); 

distributionHistories = nan(500,length(Part5Distribution.distribution().xEval)); 
evalHistories = nan(500, length(Part5Distribution.distribution().xEval));

figure(Name="UGMPDF Plot Over Time")
ax = gca; 

hold(ax, "on")

numItems = 500;
for idx = 1:numItems
    [meanVal, cov] = Part5Distribution.getCondProps(); 

    meanHistory(1, idx) = meanVal; 
    covHistory(1, idx)  = cov; 

    Zk = nonLinearState2Measurement(randomTruthSignals(1,idx));

    Part5Distribution.predictPDF_Unscented(@nonLinearStatePropNoNoise, Pww);
    Part5Distribution.updatePDF_Unscented(@nonLinearState2Measurement, Pvv,Zk); 

    [meanVal, cov] = Part5Distribution.getCondProps(); 

    meanHistory(2, idx) = meanVal; 
    covHistory(2, idx)  = cov; 


    if any(idx==1:25:500)
        dist = Part5Distribution.distribution();  
        distributionHistories(idx, :) = dist.distributions(1,:);
        evalHistories(idx, :) = dist.xEval;

        if (idx == 1)
            plot3(ax, dist.xEval(1,:), idx*ones(size(dist.xEval(1,:))), dist.distributions(1,:), "Color",'b',DisplayName="GM PDF")
        else
            plot3(ax, dist.xEval(1,:), idx*ones(size(dist.xEval(1,:))), dist.distributions(1,:), "Color",'b',DisplayName="GM PDF",HandleVisibility="off")
        end

    end
end

grid(ax, "on")

meanHistory     = reshape(meanHistory, 1,1000); 
covHistory      = reshape(covHistory, 1,1000);
arbitraryTime   = reshape([1:500; 1:500], 1, 1000); 

plot3(ax, meanHistory, arbitraryTime, zeros(1,numItems*2), DisplayName="UGMF Mean History", Color="#808080", LineWidth=1.5, LineStyle="-")

plot3(ax, meanHistory + sqrt(covHistory)*3, arbitraryTime, zeros(1,numItems*2), DisplayName="$\pm3\sigma$", Color="red", LineWidth=1.5)
plot3(ax, meanHistory - sqrt(covHistory)*3, arbitraryTime, zeros(1,numItems*2), HandleVisibility="off", Color="red", LineWidth=1.5)

plot3(ax, reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), arbitraryTime, zeros(1,numItems*2), DisplayName="True State", Color="black", LineWidth=1.5)
view(ax, 25,60)
legend(ax,Interpreter= 'LaTeX',FontSize=11)

xlabel("$x$",Interpreter= 'LaTeX',FontSize=11)
ylabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
zlabel("$p(x)$",Interpreter= 'LaTeX',FontSize=11)

figure2pdf("UGMF vs True State pdf with updates.pdf")
%%
figure(Name="Unscented GMF")
tiledlayout(2,1);

nexttile
plot(arbitraryTime, reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="True State", Color="black", LineWidth=2)
hold on
plot(arbitraryTime,meanHistory, DisplayName="UGMF Mean History", Color="#808080", LineWidth=2)
ylabel("Value",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")

nexttile
plot(arbitraryTime, meanHistory-reshape([randomTruthSignals(1,:); randomTruthSignals(1,:)], 1,numItems*2), DisplayName="UGMF Mean Error",Color="#808080", LineWidth=2, LineStyle="-")
hold on
plot(arbitraryTime,sqrt(covHistory)*3, DisplayName="$\pm3\sigma$", Color="red", LineWidth=2)
plot(arbitraryTime,- sqrt(covHistory)*3,HandleVisibility="off", Color="red", LineWidth=2)

xlabel("Time Step",Interpreter= 'LaTeX',FontSize=11)
ylabel("Error",Interpreter= 'LaTeX',FontSize=11)
legend(Interpreter= 'LaTeX',FontSize=11)
grid("on")
% 
figure2pdf("UGMF vs True State with updates.pdf")

%%
UGMF_rmse = rms(meanHistory-randomTruthSignals(1,:))


%%

% figure2pdf("UGMF vs True State with updates.pdf")

