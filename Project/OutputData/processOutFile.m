function processOutFile()
    out = evalin("base", "out");

    eulerAngles               = out.EuAngles.Data;
    quatAttitude              = out.QuatAttitude.Data; 
    
    inertialVelocity          = out.Ve.Data;
    inertialPosition          = out.Xe.Data; 
    inertialAcceleration      = out.inertialAcceleration.Data; 

    MeasuredBodyAccelerations     = squeeze(out.Ab_est.Data)'; 
    MeasuredBodyAngularRates      = squeeze(out.Wb_est.Data)';
    MeasuredBodyMagMeasurements   = squeeze(out.bodyFrameMag.Data)';     
    MeasuredVelocity              = squeeze(out.ve_est.Data)';
    MeasuredPosition              = squeeze(out.xe_est.Data)';
    
    RealBodyAccelerations     = squeeze(out.Ab_real.Data)';
    RealBodyAngularRates      = out.Wb_real.Data;
    
    BodyFrameGravity          = squeeze(out.bodyFrameGravity.Data)'; 

    % STATE_VECTOR 
    %  [1:4]  quaternion attitude
    %  [5:7]  velocity 
    %  [8:10] position
    % [11:13] gyro  bias
    % [14:16] accel bias

    TrueStateHistory          = [quatAttitude'; inertialVelocity'; inertialPosition']'; 
    
    time              = out.tout;
    
    saveName = datedFilename('Log', '.mat', ...
        'C:\Users\emmat\Desktop\ACADEMIC\YEAR5\AERO626\Project\OutputData\',"MMdd_HHmmss");

    disp("File saved at " + saveName);
    
    save(saveName,  'MeasuredBodyAccelerations', 'MeasuredBodyAngularRates', ...
                    'RealBodyAccelerations', 'RealBodyAngularRates', ...
                    'eulerAngles', 'inertialVelocity', 'inertialPosition', ...
                    "quatAttitude", "time", "MeasuredBodyMagMeasurements", ...
                    "TrueStateHistory", "inertialAcceleration", "BodyFrameGravity", ...
                    "MeasuredVelocity", "MeasuredPosition");
end