clc 
clear

% load("C:\Users\emmat\Desktop\Year5\AERO626\Project\OutputData\Log_1126_141835_SpinningYaw.mat")
% load("C:\Users\emmat\Desktop\Year5\AERO626\Project\OutputData\Log_1126_132636.mat") % Flat on Earth 
load("C:\Users\emmat\Desktop\Year5\AERO626\Project\OutputData\Log_1126_235658.mat") % Regular Flight

% STATE_VECTOR 
%  [1:4]  quaternion attitude
%  [5:7]  velocity 
%  [8:10]  position
% [11:13] gyro  bias
% [14:16] accel bias

initQuat = eul2quat(deg2rad([0 90 0]), "ZYX");

initalestimate = [initQuat 0 0 0 0 0 1 0 0 0 0 0 0];
intialestimateCovariance=eye(15).*.01;

gyro_cov    = diag([.05 .05 .05]);
accel_cov   = diag([.05 .05 .05]);
dt          = .0025; 

filter = MEKF_CLASS(initalestimate', intialestimateCovariance, gyro_cov, gyro_cov, gyro_cov, gyro_cov, .02);

[state, covariance] = filter.getFilterState();

stateHistory        = state; 
covarianceHistory   = covariance;
delf                = []; 

% trueState  = eul2quat(deg2rad([0 85 0]), "XYZ")
% firstState = state(1:4)'

firstState = rad2deg(quat2eul(state(1:4)'))




% for i=1:length(time)
numSamplesTested = 400 * 4; 
for i=1:numSamplesTested
% for i=1:100
    % angular rates are measured as x y z rates in the body frame
    % filter.propagate(MeasuredBodyAngularRates(i,:)', MeasuredBodyAccelerations(i,:)', dt); 

    % Conditional 
    if norm(MeasuredBodyAccelerations(i,:)) < 11 && norm(MeasuredBodyAccelerations(i,:)) < 9
        del = filter.updateWithGravity(MeasuredBodyAccelerations(i,:)'); 
    % elseif norm(MeasuredBodyAccelerations(i,:)) > 11
    %     disp("Norm of acceleration too high" + i)
    % 
    % elseif norm(MeasuredBodyAccelerations(i,:)) < 9
    %     disp("Norm of acceleration too low" + i)
    end
    % 
    % if mod(i, 4) == 0
    %     filter.updateWithMagnetometer(MeasuredBodyMagMeasurements(int16(i/4),:)');
    % end
    
    delf = [delf del];

    [state, covariance] = filter.getFilterState();

    % state(1:4)'

    rad2deg(quat2eul(state(1:4)'));
    
    stateHistory        = [stateHistory state]; 
    covarianceHistory   = [covariance];
end


finalError = eulerAngles(1600, :) - rad2deg(quat2eul(state(1:4)'))

figure("Name","Attitude Over Time")
estimatedTime = time(1:numSamplesTested);
EstimatedAttitude = rad2deg(quat2eul(stateHistory(1:4,1:numSamplesTested)',"ZYX"));
% RealAttitude = rad2deg(quat2eul(quatAttitude(1:numSamplesTested, :),"ZYX"));
RealAttitude = rad2deg(quat2eul(quatAttitude(1:numSamplesTested, :),"ZYX"));
tiledlayout(3,1); 
nexttile
plot(estimatedTime,EstimatedAttitude(:,1), DisplayName="Estimated Yaw");
hold on 
plot(estimatedTime,RealAttitude(:,1))
legend()
ylabel("Yaw")


nexttile
plot(estimatedTime,EstimatedAttitude(:,2), DisplayName="Estimated Pitch");
hold on 
plot(estimatedTime,RealAttitude(:,2))
legend()
ylabel("Pitch")

nexttile
plot(estimatedTime,EstimatedAttitude(:,3), DisplayName="Estimated Roll");
hold on 
plot(estimatedTime,RealAttitude(:,3))
legend()
ylabel("Roll")


figure(Name="Delf")
tiledlayout(4,1)
nexttile 
plot(estimatedTime, delf(1, : ))


nexttile 
plot(estimatedTime, delf(2, : ))


nexttile 
plot(estimatedTime, delf(3, : ))

nexttile 
plot(vecnorm(delf, 1, 1))

% figure("Name","Position Over Time")
% tiledlayout(3,1); 
% nexttile
% plot(stateHistory(8,:))
% 
% nexttile
% plot(stateHistory(9,:))
% 
% nexttile
% plot(stateHistory(10,:))

% %%


% 
% %%
% figure(Name = "Accelerations")
% 
% tiledlayout(3,1); 
% nexttile
% plot(MeasuredBodyAccelerations(1:400*4,1));
% ylabel("X Accel")
% yline(0)
% 
% nexttile
% plot(MeasuredBodyAccelerations(1:400*4,2));
% ylabel("Y Accel")
% yline(0)
% 
% nexttile
% plot(MeasuredBodyAccelerations(1:400*4,3));
% ylabel("Z Accel")
% yline(9.81)