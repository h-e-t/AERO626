clc 
clear

rng(626); 

% Problem 2 

%% Part a

Pvv = 0.02; 
Pww = 0.01^2; 

randomTruth = zeros(1,500); 
measurements = zeros(1,500); 

randomTruth(1) = 1.5; 

for idx = 2:length(randomTruth)
    randomTruth(idx)    = randomTruth(idx-1) - 0.01*sin(randomTruth(idx-1)) ...
    + randn()*chol(Pww)';
end

measurements(1) = .5*sin(2*randomTruth(1)) + randn()*chol(Pvv)';
for idx = 2:length(randomTruth)
    measurements(idx)   = 0.5*sin(2*randomTruth(idx)) + randn()*chol(Pvv)';
end


figure(Name="Truth vs Measurements")

plot(randomTruth,"black", DisplayName="Random Truth")
hold on 
plot(measurements, "rx", DisplayName="Measurements")

legend()
ylabel("Value")
xlabel("Step")

%% Part b

function nextStep = prop(x)
    Pww = 0.01^2; 
    nextStep = x - 0.01 * sin(x) + randn()*chol(Pww)'; 
end

function nextStep = propNoNoise(x)
    nextStep = x - 0.01 * sin(x); 
end

function meanMeas = state2measurement(x)
    meanMeas = .5 * sin(2*x); 
end

function F_X = Fx(x)
    F_X = 1 - .01*cos(x); 
end

Fw = 1; 

function H_X = Hx(x)
    H_X = cos(2*x); 
end

Hv = 1; 
%%
ekf_innovation_history  = [];
ekf_error_history       = [];
ekf_estimate_history    = [];
ekf_cov_history         = [];

doubleStep = []; 
regularStep = 1:500; 

mx0 = 1.5; 
pxx0 =.15^2; 

xk = mx0; 
pxxk = pxx0; 


for idx = 1:length(measurements)

    
    % Propagate
    mxkm    = prop(xk);
    pxxkm   =  Fx(xk) * pxxk * Fx(xk)' + Fw * Pww * Fw'; 
    
    % Transform state to measurement and measure
    zk = measurements(idx);
    mzkm = state2measurement(mxkm); 
    
    % Cross Covariance 
    pxzkm = pxxkm * Hx(mxkm)'; 
    
    % Innovation Covariance 
    pzzkm = Hx(mxkm) * pxxkm * Hx(mxkm)' + Hv*Pvv*Hv'; 

    % Gain 
    Kk = pxzkm/(pzzkm); 

    % Mean Update 
    mxkp = mxkm + Kk * (zk - mzkm); 

    % Covariance Update
    pxxkp = pxxkm - pxzkm*Kk' - Kk*pxzkm' + Kk * pzzkm * Kk'; 

    xk = mxkp; 
    pxxk = pxxkp; 
    

    doubleStep                  = [doubleStep idx idx];
    ekf_innovation_history      = [ekf_innovation_history pzzkm]; 
    ekf_error_history           = [ekf_error_history, 
                                        randomTruth(idx) - mxkm, 
                                        randomTruth(idx) - xk]; 
    ekf_estimate_history        = [ekf_estimate_history, mxkm, xk]; 
    ekf_cov_history             = [ekf_cov_history, pxxkm, pxxk]; 
end


figure(Name="EKF Estimation Plots")

tiledlayout(3,1) 
nexttile
plot(doubleStep, ekf_estimate_history, DisplayName="EKF Estimate History", ...
    Color="red")
hold on
plot(randomTruth, DisplayName="Truth Signal", Color="black")
legend()

nexttile 
plot(doubleStep, ekf_error_history, "black", DisplayName="Error History")
hold on

plot(doubleStep, 3*sqrt(ekf_cov_history), 'red', ...
    DisplayName="\pm3\sigma Interval")

plot(doubleStep, -3*sqrt(ekf_cov_history), 'red', ...
    DisplayName="\pm3\sigma Interval", HandleVisibility="off")
legend()

nexttile 
plot(ekf_innovation_history, 'black', DisplayName="Innovation History")
legend()




%% Part c


ukf_innovation_history_1  = []; 
ukf_cov_history_1         = [];
ukf_error_history_1       = [];
ukf_estimate_history_1    = [];

ukf_innovation_history_2  = []; 
ukf_cov_history_2         = []; 
ukf_error_history_2       = []; 
ukf_estimate_history_2    = [];

mx0 = 1.5; 
pxx0 =.15^2; 

for run = 1:2
    xk = mx0; 
    pxxk = pxx0; 
    
    switch(run)
        case 1
            alpha = 1;
            beta  = 0;
            kappa = 2;
        case 2 
            alpha = .5;
            beta  = 2;
            kappa = 2;
    end

    n       = 1;

    lambda  = alpha^2 * (n + kappa) - n;
    
    npl     = n + lambda;
        
    mean_weights    = [lambda/npl, 1/(2*npl),1/(2*npl)];
    cova_weights    = [lambda/npl + 1-alpha^2+beta, 1/(2*npl),1/(2*npl)];


    for idx = 1:length(measurements)
        % Generate Sigma Points (k-1)
        sxxkm1 = chol(pxxk, "lower");
        
        propagationSigmaPoints = [xk, xk + sqrt(npl)*sxxkm1, xk - sqrt(npl)*sxxkm1];
        
        % Propagated Sigma Points 
        chi = propNoNoise(propagationSigmaPoints);
    
        % Propagate Mean and Covariance
        mxkm    = sum(mean_weights.*chi);
    
        pxxkm = Pww;
        for jdx = 1:length(chi)
            pxxkm = pxxkm + cova_weights(jdx) * (chi(jdx) - mxkm)*(chi(jdx) - mxkm)';
        end
    
    
        % Update with measurement
        sxxkm = chol(pxxkm, "lower");
        measurementSigmaPoints     = [mxkm, mxkm + sqrt(npl)*sxxkm, mxkm - sqrt(npl)*sxxkm];
        
        
        % Transform state to measurement and measure
        zk = measurements(idx);
    
        zeta = state2measurement(measurementSigmaPoints);
        mzkm = sum(zeta.*mean_weights); 
        
        % Cross Covariance 
        pxzkm = 0;
        for jdx = 1:length(chi)
            pxzkm = pxzkm + cova_weights(jdx) * (chi(jdx) - mxkm)*(zeta(jdx) - mzkm)'; 
        end
        
        % Innovation Covariance 
        pzzkm = Pvv;
        for jdx = 1:length(chi)
            pzzkm = pzzkm + cova_weights(jdx) * (zeta(jdx) - mzkm)*(zeta(jdx) - mzkm)'; 
        end
    
        % Gain 
        Kk = pxzkm/pzzkm; 
    
        % Mean Update 
        mxkp = mxkm + Kk * (zk - mzkm); 
    
        % Covariance Update
        pxxkp = pxxkm - pxzkm*Kk' - Kk*pxzkm' + Kk * pzzkm * Kk'; 
    

        % Setup next step
        xk = mxkp; 
        pxxk = pxxkp; 
        
        switch(run)
            case 1
                ukf_innovation_history_1 = [ukf_innovation_history_1 pzzkm]; 
                ukf_error_history_1      = [ukf_error_history_1, randomTruth(idx) - mxkm, randomTruth(idx) - xk]; 

                ukf_estimate_history_1   = [ukf_estimate_history_1, mxkm, xk]; 
                ukf_cov_history_1        = [ukf_cov_history_1, pxxkm, pxxk];                
            case 2
                ukf_innovation_history_2 = [ukf_innovation_history_2 pzzkm]; 
                ukf_error_history_2      = [ukf_error_history_2, randomTruth(idx) - mxkm, randomTruth(idx) - xk]; 

                ukf_estimate_history_2   = [ukf_estimate_history_2, mxkm, xk]; 
                ukf_cov_history_2        = [ukf_cov_history_2, pxxkm, pxxk];              
        end
    end

end


figure(Name="UKF Estimation Plots")

tiledlayout(3,1) 
nexttile
plot(randomTruth, DisplayName="Truth Signal", Color="black")
hold on
plot(doubleStep, ukf_estimate_history_1, DisplayName="UKF Estimate History (i)", Color="red")
plot(doubleStep, ukf_estimate_history_2, DisplayName="UKF Estimate History (ii)", Color="#808080")
ylabel("Value")

legend()

nexttile 
plot(doubleStep, ukf_error_history_1, "black", DisplayName="Error History (i)")
hold on
plot(doubleStep, 3*sqrt(ukf_cov_history_1), 'red', DisplayName="\pm3\sigma Interval (i)")
plot(doubleStep, ukf_error_history_2, Color = "#808080", DisplayName="Error History (ii)")

plot(doubleStep, -3*sqrt(ukf_cov_history_1), 'red', HandleVisibility="off")

plot(doubleStep, 3*sqrt(ukf_cov_history_2), Color ='#FFA500', DisplayName="\pm3\sigma Interval (ii)")
plot(doubleStep, -3*sqrt(ukf_cov_history_2), Color ='#FFA500', HandleVisibility="off")

ylabel("Value")
xlabel("Step")
legend()

nexttile 
plot(ukf_innovation_history_1, 'black', DisplayName="Innovation History (i)")
hold on
plot(ukf_innovation_history_2, color = '#808080', DisplayName="Innovation History (ii)")
legend()

%% Part d combined with part e 

average_ekf_filter_error    = zeros(1,500);
average_ekf_RMS_error       = zeros(1,500); 
average_ekf_filter_variance = zeros(1,500);

average_ukf1_filter_error    = zeros(1,500);
average_ukf1_RMS_error       = zeros(1,500); 
average_ukf1_filter_variance = zeros(1,500);

average_ukf2_filter_error    = zeros(1,500);
average_ukf2_RMS_error       = zeros(1,500); 
average_ukf2_filter_variance = zeros(1,500);

for run = 1:500 
    
    randomTruth = zeros(1,500); 
    measurements = zeros(1,500); 
    
    randomTruth(1) = 1.5; 
    measurements(1) = .5*sin(2*randomTruth(1)) + randn()*chol(Pvv)';
    
    for idx = 2:length(randomTruth)
        randomTruth(idx)    = randomTruth(idx-1) - 0.01*sin(randomTruth(idx-1)) + randn()*chol(Pww)';
    end

    for idx = 2:length(randomTruth)
        measurements(idx)   = 0.5*sin(2*randomTruth(idx)) + randn()*chol(Pvv)';
    end

    ekf_error_history        = zeros(1,500); 
    ekf_cov_history          = zeros(1,500); 
    

    ukf1_error_history       = zeros(1,500); 
    ukf1_cov_history         = zeros(1,500); 

    ukf2_error_history       = zeros(1,500); 
    ukf2_cov_history         = zeros(1,500); 

    %% EKF TESTING

    
    mx0 = 1.5; 
    pxx0 =.15^2; 
    
    xk = mx0; 
    pxxk = pxx0; 

    for idx = 1:length(measurements)
        % Propagate
        mxkm    = prop(xk);
        pxxkm   =  Fx(xk) * pxxk * Fx(xk)' + Fw * Pww * Fw'; 
        
        % Transform state to measurement and measure
        zk = measurements(idx);
        mzkm = state2measurement(mxkm); 
        
        % Cross Covariance 
        pxzkm = pxxkm * Hx(mxkm)'; 
        
        % Innovation Covariance 
        pzzkm = Hx(mxkm) * pxxkm * Hx(mxkm)' + Hv*Pvv*Hv'; 
    
        % Gain 
        Kk = pxzkm/(pzzkm); 
    
        % Mean Update 
        mxkp = mxkm + Kk * (zk - mzkm); 
    
        % Covariance Update
        pxxkp = pxxkm - pxzkm*Kk' - Kk*pxzkm' + Kk * pzzkm * Kk'; 
    
        xk = mxkp; 
        pxxk = pxxkp; 
        
        ekf_error_history(idx) = randomTruth(idx) - xk; 
        ekf_cov_history(idx) = pxxk; 
    end

    %% UKF 1
    xk = mx0; 
    pxxk = pxx0; 

    alpha = 1;
    beta  = 0;
    kappa = 2;

    n       = 1;
    lambda  = alpha^2 * (n + kappa) - n; 
    npl     = n + lambda;
   
    mean_weights    = [lambda/npl, 1/(2*npl),1/(2*npl)];
    cova_weights    = [lambda/npl + 1-alpha^2+beta, 1/(2*npl),1/(2*npl)];


     for idx = 1:length(measurements)
        % Generate Sigma Points (k-1)
        sxxkm1 = chol(pxxk, "lower");
        
        propagationSigmaPoints     = [xk, xk + sqrt(npl)*sxxkm1, xk - sqrt(npl)*sxxkm1];
        
        % Propagated Sigma Points 
        chi = propNoNoise(propagationSigmaPoints);
    
        % Propagate Mean and Covariance
        mxkm    = sum(mean_weights.*chi);
    
        pxxkm = Pww;
        for jdx = 1:length(chi)
            pxxkm = pxxkm + cova_weights(jdx) * (chi(jdx) - mxkm)*(chi(jdx) - mxkm)';
        end
    
    
        % Update with measurement
        sxxkm = chol(pxxkm, "lower");
        measurementSigmaPoints     = [mxkm, mxkm + sqrt(npl)*sxxkm, mxkm - sqrt(npl)*sxxkm];
        
        
        % Transform state to measurement and measure
        zk = measurements(idx);
    
        zeta = state2measurement(measurementSigmaPoints);
        mzkm = sum(zeta.*mean_weights); 
        
        % Cross Covariance 
        pxzkm = 0;
        for jdx = 1:length(chi)
            pxzkm = pxzkm + cova_weights(jdx) * (chi(jdx) - mxkm)*(zeta(jdx) - mzkm)'; 
        end
        
        % Innovation Covariance 
        pzzkm = Pvv;
        for jdx = 1:length(chi)
            pzzkm = pzzkm + cova_weights(jdx) * (zeta(jdx) - mzkm)*(zeta(jdx) - mzkm)'; 
        end
    
        % Gain 
        Kk = pxzkm/pzzkm; 
    
        % Mean Update 
        mxkp = mxkm + Kk * (zk - mzkm); 
    
        % Covariance Update
        pxxkp = pxxkm - pxzkm*Kk' - Kk*pxzkm' + Kk * pzzkm * Kk'; 
    

        % Setup next step
        xk = mxkp; 
        pxxk = pxxkp; 

        ukf1_error_history(idx) = randomTruth(idx) - xk; 
        ukf1_cov_history(idx)   = pxxk; 
     end



    %% UKF RUN 2
    xk = mx0; 
    pxxk = pxx0;

    alpha   = .5;
    beta    = 2;
    kappa   = 2;

    n       = 1;
    lambda  = alpha^2 * (n + kappa) - n; 
    npl     = n + lambda;
       
    mean_weights    = [lambda/npl, 1/(2*npl),1/(2*npl)];
    cova_weights    = [lambda/npl + 1-alpha^2+beta, 1/(2*npl),1/(2*npl)];

    for idx = 1:length(measurements)
        % Generate Sigma Points (k-1)
        sxxkm1 = chol(pxxk, "lower");
        
        propagationSigmaPoints     = [xk, xk + sqrt(npl)*sxxkm1, xk - sqrt(npl)*sxxkm1];
        
        % Propagated Sigma Points 
        chi = propNoNoise(propagationSigmaPoints);
    
        % Propagate Mean and Covariance
        mxkm    = sum(mean_weights.*chi);
    
        pxxkm = Pww;
        for jdx = 1:length(chi)
            pxxkm = pxxkm + cova_weights(jdx) * (chi(jdx) - mxkm)*(chi(jdx) - mxkm)';
        end
    
    
        % Update with measurement
        sxxkm = chol(pxxkm, "lower");
        measurementSigmaPoints     = [mxkm, mxkm + sqrt(npl)*sxxkm, mxkm - sqrt(npl)*sxxkm];
        
        
        % Transform state to measurement and measure
        zk = measurements(idx);
    
        zeta = state2measurement(measurementSigmaPoints);
        mzkm = sum(zeta.*mean_weights); 
        
        % Cross Covariance 
        pxzkm = 0;
        for jdx = 1:length(chi)
            pxzkm = pxzkm + cova_weights(jdx) * (chi(jdx) - mxkm)*(zeta(jdx) - mzkm)'; 
        end
        
        % Innovation Covariance 
        pzzkm = Pvv;
        for jdx = 1:length(chi)
            pzzkm = pzzkm + cova_weights(jdx) * (zeta(jdx) - mzkm)*(zeta(jdx) - mzkm)'; 
        end
    
        % Gain 
        Kk = pxzkm/pzzkm; 
    
        % Mean Update 
        mxkp = mxkm + Kk * (zk - mzkm); 
    
        % Covariance Update
        pxxkp = pxxkm - pxzkm*Kk' - Kk*pxzkm' + Kk * pzzkm * Kk'; 
    

        % Setup next step
        xk = mxkp; 
        pxxk = pxxkp; 

        ukf2_error_history(idx) = randomTruth(idx) - xk; 
        ukf2_cov_history(idx)   = pxxk; 
     end



    average_ekf_filter_error(run)    = mean(ekf_error_history); 
    average_ekf_filter_variance(run) = mean(ekf_cov_history); 
    average_ekf_RMS_error(run)       =  rms(ekf_error_history); 

    average_ukf1_filter_error(run)    = mean(ukf1_error_history); 
    average_ukf1_filter_variance(run) = mean(ukf1_cov_history); 
    average_ukf1_RMS_error(run)       =  rms(ukf1_error_history); 

    average_ukf2_filter_error(run)    = mean(ukf2_error_history); 
    average_ukf2_filter_variance(run) = mean(ukf2_cov_history); 
    average_ukf2_RMS_error(run)       =  rms(ukf2_error_history); 

end

total_average_error_ekf     = mean(average_ekf_filter_error);
total_average_variance_ekf  = mean(average_ekf_filter_variance);
total_average_rmse_ekf      = mean(average_ekf_RMS_error);

total_average_error_ukf1     = mean(average_ukf1_filter_error);
total_average_variance_ukf1  = mean(average_ukf1_filter_variance);
total_average_rmse_ukf1      = mean(average_ukf1_RMS_error);

total_average_error_ukf2     = mean(average_ukf2_filter_error);
total_average_variance_ukf2  = mean(average_ukf2_filter_variance);
total_average_rmse_ukf2      = mean(average_ukf2_RMS_error);


disp(" ")
disp(" EKF mean error = " + total_average_error_ekf)
disp("UKF1 mean error = " + total_average_error_ukf1)
disp("UKF2 mean error = " + total_average_error_ukf2)
disp(" ")
disp(" EKF mean rmse = " + total_average_rmse_ekf)
disp("UKF1 mean rmse = " + total_average_rmse_ukf1)
disp("UKF2 mean rmse = " + total_average_rmse_ukf2)
disp(" ")
disp(" EKF mean Variance = " + total_average_variance_ekf)
disp("UKF1 mean Variance = " + total_average_variance_ukf1)
disp("UKF2 mean Variance = " + total_average_variance_ukf2)
% 
%  EKF mean error    = 4.6911e-05
% UKF1 mean error    = -7.0474e-05
% UKF2 mean error    = -7.1727e-05
% 
%  EKF mean rmse     = 0.051319
% UKF1 mean rmse     = 0.040561
% UKF2 mean rmse     = 0.040553
% 
%  EKF mean Variance = 0.0017325
% UKF1 mean Variance = 0.0017872
% UKF2 mean Variance = 0.0017803
% 
% The performance between the two sets of UKF parameters is nearly identical 
% with differences in mean error and rms error only apparent at the 5th 
% decimal place. 
% 
% The EKF slightly outperforms the UKF formulations in terms of mean error
% but both UKFs exhibit an improved RMSE compared to the EKF, with the second\
% set of parameters showing marginal gains. 
% 
% Final state variance slightly better in the case of the EKF with the
% second set of UKF parameters edging out the first. 
% 
% Overall, all filters illustrate practically unbiased behavior and I think
% i would choose the UKF with the second set of parameters for filtering
% measurements in this system. 


%%

figure(Name="Average Error and pm3sigma")
plot(average_ekf_filter_error, "black",DisplayName="EKF Error History")
hold on
plot(average_ukf1_filter_error,DisplayName="UKF1 Error History")
plot(average_ukf2_filter_error,DisplayName="UKF2 Error History")

plot( 3.*sqrt(average_ekf_filter_variance), "red",DisplayName="EKF \pm3\sigma ")
plot(-3.*sqrt(average_ekf_filter_variance), "red","HandleVisibility","off")

plot( 3.*sqrt(average_ukf1_filter_variance), Color="#117733",DisplayName="UKF1 \pm3\sigma ")
plot(-3.*sqrt(average_ukf1_filter_variance), Color="#117733",HandleVisibility="off")

plot( 3.*sqrt(average_ukf2_filter_variance), "red",DisplayName="UKF2\pm3\sigma ")
plot(-3.*sqrt(average_ukf2_filter_variance), "red","HandleVisibility","off")

legend()





