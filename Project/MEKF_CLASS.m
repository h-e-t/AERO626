classdef MEKF_CLASS < handle

    properties
        estimate
        estimate_covariance

        gyro_cov
        gyro_bias
        gyro_bias_cov

        accel_cov
        accel_bias
        accel_bias_cov

        barometer_variance
        barometer_bias

        previousVel
        previousAccel
        
        mag_cov
    end

    methods
        % STATE_VECTOR 
        %  [1:4]  quaternion attitude
        %  [5:7]  velocity 
        %  [8:10] position
        % [11:13] gyro  bias
        % [14:16] accel bias
        function obj = MEKF_CLASS(init_estimate, estimate_covariance, ...
                                  gyro_cov,      gyro_bias_cov,  ...
                                  accel_cov,     accel_bias_cov, baro_cov)

            obj.estimate            = init_estimate;
            obj.estimate_covariance = estimate_covariance * eye(15); 

            obj.gyro_cov            = gyro_cov;
            obj.gyro_bias           = zeros(3, 1);
            obj.gyro_bias_cov       = gyro_bias_cov; 

            obj.accel_cov           = accel_cov; 
            obj.accel_bias          = zeros(3,1); 
            obj.accel_bias_cov      = accel_bias_cov; 

            obj.mag_cov             = 1*eye(3,3); 
            
            obj.barometer_variance  = baro_cov; 
        end

        function Qd = processCovariance(obj, dt)
            Z33 = zeros(3,3); 
            
            sw = obj.gyro_cov^2; 
            sf = obj.accel_cov^2; 

            sbw = obj.gyro_bias_cov^2; 
            sbf = obj.accel_bias_cov ^ 2; 
            
            Qd = ... 
            [(sbw*dt^3)/3 + sw*dt,                         Z33,                        Z33, -(dt^2*sbw)/2,           Z33;
                              Z33,       (sbf*dt^3)/3 + sf*dt,  (sbf*dt^4)/8 + (sf*dt^2)/2,           Z33, -(dt^2*sbf)/2;
                              Z33, (sbf*dt^4)/8 + (sf*dt^2)/2, (sbf*dt^5)/20 + (sf*dt^3)/3,           Z33, -(dt^3*sbf)/6;
                   -(dt^2*sbw)/2,                         Z33,                         Z33,        dt*sbw,           Z33;
                              Z33,              -(dt^2*sbf)/2,               -(dt^3*sbf)/6,           Z33,        dt*sbf];
        end

        function obj = propagate(obj, gyro_meas, acc_meas, dt)
            gyro_meas   = gyro_meas - obj.gyro_bias;
            accel_meas  = acc_meas  - obj.accel_bias;

            assert(any(size(gyro_meas) ==[3,1]))
            assert(any(size(accel_meas)==[3,1]))
            
            %TODO: use previous values to average the integration for new
            % state estimates 
            

            % Propagate quaternion attiude using measurement
            quatAttitude = quaternion(obj.estimate(1:4)') * quaternion(quatexp(.5*reshape([0;gyro_meas.*dt],1,4)));

            % update the velocity and position estimates using quaternion
            %              (body2intertialrotation * measured acceleration + gravity) * dt + previous velocity 
            vkm          = (quat2rotm(quatAttitude)*accel_meas + [0 0 9.81]').*dt + obj.estimate(5:7); 
            rkm          = vkm*dt + obj.estimate(8:10); 
            
            % Covariance Propagation 
            I3 = eye(3); 
            Z3 = zeros(3,3); 

            %TODO: check about specific force vs acceleration in F
            F  = [-skewSym(gyro_meas)                          Z3 Z3 -I3 Z3                      ; 
                  -quat2rotm(quatAttitude)*skewSym(accel_meas) Z3 Z3  Z3 -quat2rotm(quatAttitude);
                  Z3                                           I3 Z3  Z3 Z3;
                  Z3                                           Z3 Z3  Z3 Z3;
                  Z3                                           Z3 Z3  Z3 Z3];
            
            % Calculate new covariance
            obj.estimate_covariance = F * obj.estimate_covariance * F' + obj.processCovariance(dt);

            obj.estimate(1:4)   = quatAttitude.compact; 
            obj.estimate(5:7)   = vkm; 
            obj.estimate(8:10)  = rkm;             
        end

        function delf = updateWithGravity(obj, acc_meas)
            aPrioriErrorState = zeros(15,1); % Only in attitude

            Z33 = zeros(3,3); 
            I33  = eye(3); 

            H = [skewSym(-quat2rotm(obj.estimate(1:4)')*[0;0;9.81]) , Z33, Z33, Z33, I33];

            delf = quat2rotm(obj.estimate(1:4)')'*[0;0;-9.81] - (acc_meas - obj.estimate(14:16))
            % delf = (acc_meas - obj.estimate(14:16)) - quat2rotm(obj.estimate(1:4)')'*[0;0;-9.81];
            
            Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.accel_cov);

            aPosterioriErrorState = aPrioriErrorState + Kgain*delf;
            
            % aPosterior error state covariance
            obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance; 
            
            obj.estimate(1:4) = quatmultiply(obj.estimate(1:4)', [1 aPosterioriErrorState(1:3)'./2]); 
        end

        function obj = updateWithBarometer(obj, baro_meas)
            aPrioriErrorState = zeros(15,1); % Only in attitude
            
            Z13 = zeros(1,3); 

            H = [Z13 Z13 [0 0 1] Z13 Z13]; 

            del_z  = baro_meas - obj.estimate(10);

            % TODO: Include barometer measurement variance / bias estimates 
            Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.barometer_variance);

            aPosterioriErrorState = aPrioriErrorState + Kgain*del_z

            obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance; 

            obj.estimate(2:end) = obj.estimate(2:end) + aPosterioriErrorState; 
            
        end

        function obj = updateWithMagnetometer(obj, magMeas)
            aPrioriErrorState = zeros(15,1); % Only in attitude
            
            Z33 = zeros(3,3); 
            I33  = eye(3); 

            NED_mag = [27.5550; -2.4169; -16.0849];

            H = [skewSym(-quat2rotm(obj.estimate(1:4)')*NED_mag), Z33, Z33, Z33, Z33];

            % delm = quat2rotm(obj.estimate(1:4)')'*-NED_mag - magMeas
            delm = magMeas - quat2rotm(obj.estimate(1:4)')'*-NED_mag;

            Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.mag_cov);

            aPosterioriErrorState = aPrioriErrorState + Kgain*delm;
            
            % aPosterior error state covariance
            obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance; 
            
            obj.estimate(1:4) = quatmultiply(obj.estimate(1:4)', [1 aPosterioriErrorState(1:3)'./2]); 
        end

        function [state, covariance] = getFilterState(obj)
            state       = obj.estimate; 
            covariance  = diag(obj.estimate_covariance); 
        end

    end
end
