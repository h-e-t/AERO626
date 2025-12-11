function Plotter(time, EstimatedStateHistory, EstimatedCovarianceHistory, TrueStateHistory, plotPosition, plotAttitude, plotVelocity)
    arguments 
        time
        EstimatedStateHistory 
        EstimatedCovarianceHistory 
        TrueStateHistory 
        plotPosition = 0
        plotAttitude = 0
        plotVelocity = 0
    end
    % STATE_VECTOR 
    %  [1:4]  quaternion attitude
    %  [5:7]  velocity 
    %  [8:10] position
    % [11:13] gyro  bias
    % [14:16] accel bias
    
    if plotPosition
        figure("Name","Estimated Position Over Time")
        tiledlayout(3,1); 
        nexttile
        plot(time, EstimatedStateHistory(:,8), DisplayName="Estimated X")
        hold on
        plot(time, TrueStateHistory(:,8), DisplayName="True X")
        legend()
        
        nexttile
        plot(time, EstimatedStateHistory(:,9), DisplayName="Estimated Y")
        hold on
        plot(time, TrueStateHistory(:,9), DisplayName="True Y")
        legend()
        
        nexttile
        plot(time, -EstimatedStateHistory(:,10), DisplayName="Estimated Z")
        hold on
        plot(time, -TrueStateHistory(:,10), DisplayName="True Z")
        legend()
        
    end

    if plotAttitude
        EstimatedAttitude = eulerd(quaternion(EstimatedStateHistory(:, 1:4)),"ZYX","frame");
        RealAttitude      = eulerd(quaternion(TrueStateHistory(:, 1:4)),"ZYX","frame");

        figure(Name="Attitude Over time")
        tiledlayout(3,1); 
        nexttile
        plot(time,EstimatedAttitude(:,1), DisplayName="Estimated Yaw");
        hold on 
        plot(time,RealAttitude(:,1), DisplayName="True Yaw")
        legend()
        ylabel("Yaw")
        
        
        nexttile
        plot(time,EstimatedAttitude(:,2), DisplayName="Estimated Pitch");
        hold on 
        plot(time,RealAttitude(:,2), DisplayName="True Pitch")
        legend()
        ylabel("Pitch")
        
        nexttile
        plot(time,EstimatedAttitude(:,3), DisplayName="Estimated Roll");
        hold on 
        plot(time,RealAttitude(:,3), DisplayName="True Roll")
        legend()
        ylabel("Roll")
    end

    if plotVelocity
        figure("Name","Estimated Inertial Velocity Over Time")
        tiledlayout(3,1); 
        nexttile
        plot(time, EstimatedStateHistory(:,5), DisplayName="Estimated X Vel")
        hold on
        plot(time, TrueStateHistory(:,5), DisplayName="True X Vel")
        legend()
        
        nexttile
        plot(time, EstimatedStateHistory(:,6), DisplayName="Estimated Y Vel")
        hold on
        plot(time, TrueStateHistory(:,6), DisplayName="True Y Vel")
        legend()
        
        nexttile
        plot(time, -EstimatedStateHistory(:,7), DisplayName="Estimated Z Vel")
        hold on
        plot(time, -TrueStateHistory(:,7), DisplayName="True Z Vel")
        legend()
    end
end