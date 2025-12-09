
load("OutputData\Log_1208_171059.mat")

figure(Name="Position over time")
tiledlayout(3,1)
nexttile 
plot(time, inertialPosition(:,1), "black",LineWidth=2)
hold on; 
plot(time(1:80:end), MeasuredPosition(1:80:end, 1), "red", LineWidth=2)

ylabel("X Inertial Position (m)")

nexttile 
plot(time, inertialPosition(:,2), "black",LineWidth=2)
hold on; 
plot(time(1:80:end), MeasuredPosition(1:80:end, 2), "red", LineWidth=2)
ylabel("Y Inertial Position (m)")


nexttile 
plot(time, -inertialPosition(:,3)+150, "black",LineWidth=2)
hold on; 
plot(time(1:80:end), -MeasuredPosition(1:80:end, 3)+150, "red", LineWidth=2)
ylabel("Z Inertial Position (m)")
xlabel("Time (s)")

exportgraphics(gcf,"interialPosOverTime.png",Units="inches", ...
       Width=5,Height=4,Resolution=300)

%%
figure(Name="Velocity over time")
tiledlayout(3,1)
nexttile 
plot(time, inertialVelocity(:,1), "black",LineWidth=2)
hold on; 
plot(time(1:80:end), MeasuredVelocity(1:80:end, 1), "red", LineWidth=2)
ylabel("X Inertial Velocity (m)")

nexttile 
plot(time, inertialVelocity(:,2), "black",LineWidth=2)
hold on; 
plot(time(1:80:end), MeasuredVelocity(1:80:end, 2), "red", LineWidth=2)
ylabel("Y Inertial Velocity (m)")


nexttile 
plot(time, -inertialVelocity(:,3), "black",LineWidth=2)
hold on; 
plot(time(1:80:end), -MeasuredVelocity(1:80:end, 3), "red", LineWidth=2)
ylabel("Z Inertial Velocity (m)")
xlabel("Time (s)")

exportgraphics(gcf,"interialVelOverTime.png",Units="inches", ...
       Width=5,Height=4,Resolution=300)


%%
figure("Name","MeasurementVSTrueOvertime")

tiledlayout(3,1)
nexttile 
plot(time, MeasuredBodyAngularRates(:,1), "red",LineWidth=2, DisplayName="Measured Angular Rates")
hold on
plot(time, RealBodyAngularRates(:,1), "black",LineWidth=2, DisplayName="Real Angular Rates")
legend()
ylabel("X Rate (rad/s)")

nexttile 
plot(time, MeasuredBodyAngularRates(:,2), "red",LineWidth=2)
hold on
plot(time, RealBodyAngularRates(:,2), "black",LineWidth=2)
ylabel("Y Rate (rad/s)")

nexttile 
plot(time, MeasuredBodyAngularRates(:,3), "red",LineWidth=2)
hold on
plot(time, RealBodyAngularRates(:,3), "black",LineWidth=2)
ylabel("Z Rate (rad/s)")

xlabel("Time (s)")

exportgraphics(gcf,"realvsmeasuredangularRates.png",Units="inches", ...
       Width=5,Height=4,Resolution=300)

%%
figure("Name","AccelMeasurementVSTrueOvertime")

tiledlayout(3,1)
nexttile 
plot(time, MeasuredBodyAccelerations(:,1), "red",LineWidth=2, DisplayName="Measured Acceleration")
hold on
plot(time, RealBodyAccelerations(:,1), "black",LineWidth=2, DisplayName="Real Acceleration")
legend()
ylabel("X Acceleration (m/s^2)")

nexttile 
plot(time, MeasuredBodyAccelerations(:,2), "red",LineWidth=2)
hold on
plot(time, RealBodyAccelerations(:,2), "black",LineWidth=2)
ylabel("Y Acceleration (m/s^2)")


nexttile 
plot(time, MeasuredBodyAccelerations(:,3), "red",LineWidth=2)
hold on
plot(time, RealBodyAccelerations(:,3), "black",LineWidth=2)
ylabel("Z Acceleration (m/s^2)")

xlabel("Time (s)")

exportgraphics(gcf,"realvsmeasuredaccelerations.png",Units="inches", ...
       Width=5,Height=4,Resolution=300)

