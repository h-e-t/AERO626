clc

Pww = .01^2; 

Pvv = .02; 

means = 1.5;
means = reshape(means, 1,1,[]);

covariances = .15^2; 
covariances = reshape(covariances,1,1,[]);

weights = 1;


Part5Distribution = ExtendedGaussianMixtureFilter(means, covariances, weights);

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

    Part5Distribution.predictPDF_Extended(@nonLinearStatePropNoNoise, @linearizedDynamics,Pww);
    Part5Distribution.updatePDF_Extended(@nonLinearState2Measurement, @linearizedMeasurements, Pvv,Zk); 

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

% figure2pdf("UGMF vs True State pdf with updates.pdf")
% %%
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
