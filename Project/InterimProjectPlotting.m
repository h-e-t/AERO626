figure(Name="Position over time")
tiledlayout(3,1)
nexttile 
plot(time, inertialPosition(:,1), "black",LineWidth=2)
ylabel("X Inertial Position (m)")

nexttile 
plot(time, inertialPosition(:,2), "black",LineWidth=2)
ylabel("Y Inertial Position (m)")


nexttile 
plot(time, -inertialPosition(:,3), "black",LineWidth=2)
ylabel("Z Inertial Position (m)")
xlabel("Time (s)")

exportgraphics(gcf,"interialPosOverTime.png",Units="inches", ...
       Width=5,Height=4,Resolution=300)


%%

figure("Name","MeasurementVSTrueOvertime")

tiledlayout(3,1)
nexttile 
plot(time, MeasuredBodyAngularRates(:,1), "red",LineWidth=2, DisplayName="Measured Angular Rates")
hold on
plot(time, RealBodyAngularRates(:,1), "black",LineWidth=2, DisplayName="Real Angular Rates")
legend()
ylabel("X Angular Rates (rad/s)")

nexttile 
plot(time, MeasuredBodyAngularRates(:,2), "red",LineWidth=2)
hold on
plot(time, RealBodyAngularRates(:,2), "black",LineWidth=2)
ylabel("Y Angular Rates (rad/s)")



nexttile 
plot(time, MeasuredBodyAngularRates(:,3), "red",LineWidth=2)
hold on
plot(time, RealBodyAngularRates(:,3), "black",LineWidth=2)
ylabel("Z Angular Rates (rad/s)")

xlabel("Time (s)")

exportgraphics(gcf,"realvsmeasuredangularRates.png",Units="inches", ...
       Width=5,Height=4,Resolution=300)
%%

figure("Name","AccelMeasurementVSTrueOvertime")

tiledlayout(3,1)
nexttile 
plot(time, MeasuredBodyAccelerations(:,1), "red",LineWidth=2, DisplayName="Measured Acceleration")
hold on
plot(time, RealBodyAccelerations(:,1)-BodyFrameGravity(:,1), "black",LineWidth=2, DisplayName="Real Acceleration")
legend()
ylabel("X Acceleration (m/s^2)")

nexttile 
plot(time, MeasuredBodyAccelerations(:,2), "red",LineWidth=2)
hold on
plot(time, RealBodyAccelerations(:,2)-BodyFrameGravity(:,2), "black",LineWidth=2)
ylabel("Y Acceleration (m/s^2)")


nexttile 
plot(time, MeasuredBodyAccelerations(:,3), "red",LineWidth=2)
hold on
plot(time, RealBodyAccelerations(:,3)-BodyFrameGravity(:,3), "black",LineWidth=2)
ylabel("Z Acceleration (m/s^2)")

xlabel("Time (s)")

exportgraphics(gcf,"realvsmeasuredaccelerations.png",Units="inches", ...
       Width=5,Height=4,Resolution=300)

%%

figure("Name","Euler Angles Over Time")

tiledlayout(3,1)
nexttile
angles = rad2deg(eulerAngles); 

plot(time, angles(:,1))

nexttile

plot(time, angles(:,2))

nexttile

plot(time, angles(:,3))