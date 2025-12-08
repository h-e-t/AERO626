clc 
clear

load("OutputData\Log_1206_130320.mat")
%
% STATE_VECTOR 
%  [1:4]  quaternion attitude
%  [5:7]  velocity 
%  [8:10]  position
% [11:13] gyro  bias
% [14:16] accel bias

% COVARIANCE VECTOR
%  [1:3]  quaternion covariance
%  [4:6]  velocity covariance 
%  [7:9]  position covariance
% [10:12] gyro bias covariance
% [13:15] accel bias covariance

initQuat = eul2quat(deg2rad([15 76 5]), "ZYX");

initalestimate = [initQuat 0 0 0 0 0 0 0 0 0 0 0 0];
intialestimateCovariance=eye(15).*.01;

gyro_cov    = diag([.05 .05 .05]);
gyro_bias_cov    = diag([.05 .05 .05]);

accel_cov   = diag([.05 .05 .05]);
accel_bias_cov    = diag([.05 .05 .05]);

mag_cov   = diag([.05 .05 .05]);

gps_cov   = diag([.05 .05 .05 1.17 1.17 1.17]);

baro_cov    = .5;

dt          = .0025; 


initFilter = MEKF_CLASS(initalestimate', intialestimateCovariance, ...
                    gyro_cov, gyro_bias_cov, ...
                    accel_cov, accel_bias_cov, ...
                    mag_cov, gps_cov, baro_cov);


filter = EGMF_CLASS(initFilter);


%%
[state, covariance] = filter.getFilterState();

numSamplesTested = 4300; % this completes the ascent for 85 deg 

stateHistory        = nan(numSamplesTested,16); 
covarianceHistory   = nan(numSamplesTested,15);

for i=1:numSamplesTested
    if mod(i, 200) == 0
        disp("Samples Cleared: ")
    end
    [state, covariance] = filter.getFilterState();
    
    stateHistory(i, :)        = state'; 
    covarianceHistory(i, :)   = covariance';

    BodyAccelBias = [0 0 0]';
    % angular rates are measured as x y z rates in the body frame
    filter.propagate(MeasuredBodyAngularRates(i,:)', MeasuredBodyAccelerations(i,:)' + BodyAccelBias, dt);
    
    filter.updateWithGravity(MeasuredBodyAccelerations(i,:)');

    % Magnetometer Update (Sampled lower to reflect update frequency)
    if mod(i, 4) == 0 
        filter.updateWithMagnetometer(MeasuredBodyMagMeasurements(i,:)');
    end

    % GPS Update (Sampled lower to reflect update frequency)
    if mod(i, 80) == 0  
        % Rayleight Sampling
        [pos_meas, vel_meas] = RayleighSampling(inertialVelocity(i, :)',inertialPosition(i, :)');

        filter.updateWithGPS(pos_meas, vel_meas);
    end

    % Baro Update
    if mod(i, 15) == 0
        alt = TrueStateHistory(i, 10) + baro_cov*randn(); 
        filter.updateWithBarometer(alt);
    end


end
%
%%
figure(Name = "Quaternion Estimate Covariance")
quatError = stateHistory(:, 2:4) - TrueStateHistory(1:numSamplesTested, 2:4); 
plotTiledCovariance(3,time(1:numSamplesTested), quatError, 3*sqrt(covarianceHistory(:, 1:3)),{"$q_1$ Error", "$q_2$ Error", "$q_3$ Error", "$q_4$ Error"}, "Time (s)")


figure(Name = "Velocity Estimate Covariance")
velError = stateHistory(:, 5:7) - TrueStateHistory(1:numSamplesTested, 5:7);
plotTiledCovariance(3,time(1:numSamplesTested), velError, 3*sqrt(covarianceHistory(:,4:6)),{"$v_x$ Error", "$v_y$ Error", "$v_z$ Error"}, "Time (s)")


figure(Name = "Position Estimate Covariance")
posError = stateHistory(:,8:10) - TrueStateHistory(1:numSamplesTested,8:10);
plotTiledCovariance(3,time(1:numSamplesTested), posError, 3*sqrt(covarianceHistory(:, 7:9)),{"$r_x$ Error", "$r_y$ Error", "$r_z$ Error"}, "Time (s)")

%%
Plotter(time(1:numSamplesTested), stateHistory,covarianceHistory, TrueStateHistory(1:numSamplesTested,:), 0, 1,0); 


%%
% figure(Name = "Gyro bias estimate")
% estimateError = stateHistory(:,14:16); 
% plotTiledCovariance(3,time(1:numSamplesTested), estimateError, 3*sqrt(covarianceHistory(:,10:12)),{"$b_{\omega x}$ Error", "$b_{\omega y}$ Error", "$b_{\omega z}$ Error"}, "Time (s)")

%%
% figure(Name = "Accel bias estimate")
% estimateError = stateHistory(:, 11:13); 
% plotTiledCovariance(3,time(1:numSamplesTested), estimateError, 3*sqrt(covarianceHistory(:,13:15)),{"$r_x$ Error", "$r_y$ Error", "$r_z$ Error"}, "Time (s)")

