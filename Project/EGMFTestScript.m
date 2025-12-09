clc 
clear

load("OutputData\Log_1208_211402.mat")
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

initQuat = eul2quat(deg2rad([15 75 5]), "ZYX");

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
        disp("Samples Cleared: " + i)
    end
    [state, covariance] = filter.getFilterState();
    
    stateHistory(i, :)        = state'; 
    covarianceHistory(i, :)   = covariance';

    % angular rates are measured as x y z rates in the body frame
    filter.propagate(MeasuredBodyAngularRates(i,:)', MeasuredBodyAccelerations(i,:)', dt);
    
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
figure(Name = "Velocity Estimate Covariance")

velError = stateHistory(:, 5:7) - TrueStateHistory(1:numSamplesTested, 5:7);
plotTiledCovariance(3,time(1:numSamplesTested), velError, 3*sqrt(covarianceHistory(:,4:6)),{"$v_x$ Error", "$v_y$ Error", "$v_z$ Error"}, "Time (s)")

figure2pdf("EGMF Velocity Estimate Covariance");
%%
figure(Name = "Position Estimate Covariance")
posError = stateHistory(:,8:10) - TrueStateHistory(1:numSamplesTested,8:10);
plotTiledCovariance(3,time(1:numSamplesTested), posError, 3*sqrt(covarianceHistory(:, 7:9)),{"$r_x$ Error", "$r_y$ Error", "$r_z$ Error"}, "Time (s)")

figure2pdf("EGMF Position Estimate Covariance");
% 
%%
figure(Name = "Attitude Estimate Covariance")

quatError = quatmultiply(stateHistory(:, 1:4), quatinv(TrueStateHistory(1:numSamplesTested, 1:4))); 
angles = rad2deg(quat2eul(quatError, "ZYX"));

plotTiledCovariance(3,time(1:numSamplesTested), angles, 5.*sqrt(covarianceHistory(:, 1:3)),{"Yaw Error", "Pitch Error", "Roll Error"}, "Time (s)")
figure2pdf("EGMF Attitude Estimate Covariance");
%%
Plotter(time(1:numSamplesTested), stateHistory,covarianceHistory, TrueStateHistory(1:numSamplesTested,:), 0, 1,0); 

