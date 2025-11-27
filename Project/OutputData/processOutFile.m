function processOutFile()
    out = evalin("base", "out");

    MeasuredBodyAccelerations = out.Ab_est.Data; 
    MeasuredBodyAngularRates  = out.Wb_est.Data;
    MeasuredBodyMagMeasurements   = squeeze(out.bodyFrameMag.Data)'; 
    
    BodyFrameGravity          = squeeze(out.bodyFrameGravity.Data)';
    
    
    RealBodyAccelerations     = out.Ab_real.Data;
    RealBodyAngularRates      = out.Wb_real.Data;
    
    
    eulerAngles               = out.EuAngles.Data;
    quatAttitude              = out.QuatAttitude.Data; 
    
    inertialVelocity          = out.Ve.Data;
    inertialPosition          = out.Xe.Data; 
    
    
    
    time              = out.tout;
    
    saveName = datedFilename('Log', '.mat', ...
        'C:\Users\emmat\Desktop\Year5\AERO626\Project\OutputData\',"MMdd_HHmmss");

    disp("File saved at " + saveName);
    
    save(saveName,  'MeasuredBodyAccelerations', 'MeasuredBodyAngularRates', ...
                    'RealBodyAccelerations', 'RealBodyAngularRates', ...
                    'eulerAngles', 'inertialVelocity', 'inertialPosition', ...
                    "quatAttitude", "time", "BodyFrameGravity", "BodyFrameGravity", ...
                    "MeasuredBodyMagMeasurements");
end