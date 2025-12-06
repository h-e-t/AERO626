clc 
clear

% load("OutputData\FirstTimeWithCustomAccelBlock.mat")
load("OutputData\Log_1205_173919.mat")

% TrueStateHistory = [quatAttitude inertialVelocity inertialPosition]';
MeasuredBodyAccelerations = squeeze(MeasuredBodyAccelerations)';
MeasuredBodyAngularRates =  squeeze(MeasuredBodyAngularRates)';
%%
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

initQuat = eul2quat(deg2rad([0 75 15]), "ZYX");

initalestimate = [initQuat 0 0 0 0 0 0 0 0 0 0 0 0];
intialestimateCovariance=eye(15).*.01;

gyro_cov    = diag([.05 .05 .05]);
accel_cov   = diag([.05 .05 .05]);
dt          = .0025; 

filter = MEKF_CLASS(initalestimate', intialestimateCovariance, gyro_cov, gyro_cov, gyro_cov, gyro_cov, 4);

[state, covariance] = filter.getFilterState();

stateHistory        = []; 
covarianceHistory   = [];
delf                = []; 
delm                = []; 
magdelm             = [];
delz                = [];
covDelta            = []; 
covNorm             = [];
fNorm             = [];

firstState = rad2deg(quat2eul(state(1:4)'))

numSamplesTested = 4634; % this completes the ascent for 85 deg 

for i=1:numSamplesTested
    BodyAccelBias = [.2 .2 -.2]';
    % angular rates are measured as x y z rates in the body frame
    [delta, cov] = filter.propagate(MeasuredBodyAngularRates(i,:)', MeasuredBodyAccelerations(i,:)' + BodyAccelBias, dt);
    filter.estimate_covariance = cov; 

    covDelta = [covDelta, norm(delta)];
    covNorm = [covNorm norm(cov)];
    
    % Conditional 
    if norm(MeasuredBodyAccelerations(i,:)) < 11 && norm(MeasuredBodyAccelerations(i,:)) > 8.5 && i < 2000
        del = filter.updateWithGravity(MeasuredBodyAccelerations(i,:)');
        delf = [delf del];
    end

    % Magnetometer Update
    if mod(i, 4) == 0 && i < 2000
        del = filter.updateWithMagnetometer(MeasuredBodyMagMeasurements(int16(i/4),:)');
        delm = [delm del];
        magdelm = [magdelm norm(del)];

    end
    
    if mod(i, 80) == 0
        % Rayleight Sampling
        [pos_meas, vel_meas] = RayleighSampling(inertialVelocity(i, :)',inertialPosition(i, :)');

        filter.updateWithGPS(pos_meas, vel_meas);
    end


    if mod(i, 15) == 0
        alt = TrueStateHistory(10, i); 
        del = filter.updateWithBarometer(alt);
        delz = [delz del];

    end
    
    [state, covariance] = filter.getFilterState();
    
    stateHistory        = [stateHistory state]; 
    covarianceHistory   = [covarianceHistory covariance];
end

% finalError = rad2deg(quat2eul(state(1:4)')) - rad2deg(eulerAngles(numSamplesTested, :))

% 
Plotter(time(1:numSamplesTested), stateHistory,covarianceHistory, TrueStateHistory(:,1:numSamplesTested), 1, 1, 1); 
%%
figure(Name = "Quaternion Estimate Covariance")
quatError = stateHistory(1:4,:) - TrueStateHistory(1:4,1:numSamplesTested); 
plotTiledCovariance(4,time(1:numSamplesTested), quatError, 3*sqrt(covarianceHistory(1:4,:)),{"$q_1$ Error", "$q_2$ Error", "$q_3$ Error", "$q_4$ Error"}, "Time (s)")

%%
figure(Name = "Velocity Estimate Covariance")
velError = stateHistory(5:7,:) - TrueStateHistory(5:7,1:numSamplesTested);
plotTiledCovariance(3,time(1:numSamplesTested), velError, 3*sqrt(covarianceHistory(4:6,:)),{"$v_x$ Error", "$v_y$ Error", "$v_z$ Error"}, "Time (s)")

%%
figure(Name = "Position Estimate Covariance")
posError = stateHistory(8:10,:) - TrueStateHistory(8:10,1:numSamplesTested);
plotTiledCovariance(3,time(1:numSamplesTested), posError, 3*sqrt(covarianceHistory(7:9,:)),{"$r_x$ Error", "$r_y$ Error", "$r_z$ Error"}, "Time (s)")

figure(Name = "Gyro bias estimate")
estimateError = stateHistory(14:16, :); 
plotTiledCovariance(3,time(1:numSamplesTested), estimateError, 3*sqrt(covarianceHistory(10:12,:)),{"$b_{\omega x}$ Error", "$b_{\omega y}$ Error", "$b_{\omega z}$ Error"}, "Time (s)")

figure(Name = "Accel bias estimate")
estimateError = stateHistory(11:13, :); 
plotTiledCovariance(3,time(1:numSamplesTested), estimateError, 3*sqrt(covarianceHistory(13:15,:)),{"$r_x$ Error", "$r_y$ Error", "$r_z$ Error"}, "Time (s)")


%%
figure()
tiledlayout(3,1)
for i = 1:3
    nexttile
    plot(MeasuredBodyAccelerations(1:numSamplesTested,i))
    hold on
    plot(RealBodyAccelerations(1:numSamplesTested, i)-BodyFrameGravity(1:numSamplesTested, i))
end


%%
figure()

figure()
tiledlayout(4,1)
for i = 1:4
    nexttile

    plot(stateHistory(i,:)-TrueStateHistory(i,1:numSamplesTested))
    % hold on 
    % plot(TrueStateHistory(i,1:numSamplesTested))

end

