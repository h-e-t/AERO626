clear 
clc

% load("data_HW06.mat")

F = [1 1; 
     0 1];

H = [1 0];

Pww = [0 0    ; 
       0    0.01]; 

Pvv = 1; 

means = [[-2.5 -.5]' [-1 .2]' [-.5 -.3]' [1.2 1]'];
means = reshape(means, 2,1,[]);

covariances = [diag([.2 .1]) diag([.25 .05]) diag([.2 .1]) diag([1 .3])];
covariances = reshape(covariances,2,2,[]);


chol(covariances,'lower')
%%

weights = [.2 .3 .1 .4];
 

figure; 
tiledlayout(1,2)
ax1 = nexttile(1);
hold on;
ax2 = nexttile(2);
hold on;

mixture = GaussianMixtureFilter(means, covariances, weights);
dist = mixture.distribution();
posPDF = dist.distributions(1,:); 
velPDF = dist.distributions(2,:); 

plot3(ax1, dist.xEval(1,:), zeros(size(dist.xEval(1,:))), posPDF, "Color",'b',DisplayName="GM PDF")
plot3(ax2, dist.xEval(2,:), zeros(size(dist.xEval(2,:))), velPDF, "Color",'b','HandleVisibility','off')
view(ax2, 25,60)
view(ax1, 25,60)
%%
covarianceHistory   = zeros(2,length(zk)*2+1); 
meanHistory         = zeros(2,length(zk)*2+1);
errorHistory        = zeros(2,length(zk)*2+1);

[mean, cov] = mixture.getCondProps(); 

historyIdx = 1; 

historyTimeStep = [0];
meanHistory(:,historyIdx)       = mean; 
covarianceHistory(:,historyIdx) = diag(cov); 
errorHistory(:,historyIdx)      = x0-mean; 
historyIdx = historyIdx+1; 

for idx = 1:length(zk)
% for idx = 1:6
    mixture.predictPDF(F, Pww); 
    
    % A Priori Mean and Covariance 
    [mean, cov]                         = mixture.getCondProps();
    meanHistory(:,historyIdx)           = mean; 
    errorHistory(:,historyIdx)          = mean - xk(:,idx);
    covarianceHistory(:,historyIdx)     = diag(cov); 
    
    mixture.updatePDF(H, Pvv, zk(idx)); 

    % A Posteriori Mean and Covariance 
    [mean, cov]                         = mixture.getCondProps(); 
    meanHistory(:,historyIdx+1)         = mean; 
    errorHistory(:,historyIdx+1)        = mean - xk(:,idx);
    covarianceHistory(:,historyIdx+1)   = diag(cov); 
    
    dist    = mixture.distribution();    
    posPDF = dist.distributions(1,:); 
    velPDF = dist.distributions(2,:); 

    historyTimeStep = [historyTimeStep idx idx];
    historyIdx = historyIdx+2; 
    plot3(ax1, dist.xEval(1,:), idx*ones(size(dist.xEval(1,:))), posPDF, "Color",'b','HandleVisibility','off')
    plot3(ax2, dist.xEval(2,:), idx*ones(size(dist.xEval(2,:))), velPDF, "Color",'b','HandleVisibility','off')
    
end 

plot3(ax1, xk(1,:), 1:length(zk), zeros(1,length(zk)), DisplayName="True State")
legend(ax1)
plot3(ax2, xk(2,:), 1:length(zk), zeros(1,length(zk)))


grid(ax1, "on")
xlabel(ax1, "Position", Interpreter="latex")
ylabel(ax1, "Time Step, k", Interpreter="latex")
zlabel(ax1, "$$\textit{p(x)}$$", Interpreter="latex")
view(ax1, 25,60)

grid(ax2, "on")
xlabel(ax2, "Velocity", Interpreter="latex")
ylabel(ax2, "Time Step, k", Interpreter="latex")
zlabel(ax2, "$$\textit{p(x)}$$", Interpreter="latex")
view(ax2, 25,60)


figure
tiledlayout(2,1)
nexttile
plot(historyTimeStep, errorHistory(1,:), "black", DisplayName="GMF Error History")
hold on;
plot(historyTimeStep, 3.*sqrt(covarianceHistory(1,:)), "red",DisplayName="\pm3\sigma ")
plot(historyTimeStep,-3.*sqrt(covarianceHistory(1,:)), "red","HandleVisibility","off")
ylabel("Position Error")
legend()


nexttile
plot(historyTimeStep, errorHistory(2,:), "black", DisplayName="GMF Velocity Error History")
hold on;
plot(historyTimeStep, 3.*sqrt(covarianceHistory(2,:)), "red",DisplayName="\pm3\sigma ")
plot(historyTimeStep,-3.*sqrt(covarianceHistory(2,:)), "red","HandleVisibility","off")
ylabel("Velocity Error")

xlabel("Time Step")

%%
posRMSE = rms(errorHistory(1,:))
velRMSE = rms(errorHistory(2,:))
