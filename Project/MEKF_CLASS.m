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

        propagations_since_last_update
        
        mag_cov

        gpsVariance
    end

    methods
        % STATE_VECTOR 
        %  [1:4]  quaternion attitude
        %  [5:7]  velocity 
        %  [8:10] position
        % [11:13] gyro  bias
        % [14:16] accel bias

        % ERROR STATE_VECTOR 
        %  [1:3]  quaternion attitude
        %  [4:6]  velocity 
        %  [7:9]  position
        % [10:12] gyro  bias
        % [13:15] accel bias
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

            obj.propagations_since_last_update = 0;

            obj.gpsVariance         = diag([.1 .1 .1 1.17 1.17 1.17]); 
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
                              Z33, (sbf*dt^4)/8 + (sf*dt^2)/2, (sbf*dt^5)/20 + (sf*dt^3)/3 ,           Z33, -(dt^3*sbf)/6;
                   -(dt^2*sbw)/2,                         Z33,                         Z33,        dt*sbw,           Z33;
                              Z33,              -(dt^2*sbf)/2,               -(dt^3*sbf)/6,           Z33,        dt*sbf];
        end

        function [delta, cov] = propagate(obj, gyro_meas, acc_meas, dt)
            gyro_meas   = gyro_meas - obj.gyro_bias;
            accel_meas  = acc_meas  - obj.accel_bias;

            assert(any(size(gyro_meas) ==[3,1]))
            assert(any(size(accel_meas)==[3,1]))
            
            % Propagate quaternion attiude using measurement
            if norm(obj.estimate(1:4)') ~= 1
                obj.estimate(1:4) = obj.estimate(1:4) / norm(obj.estimate(1:4)');
            end


            % update the velocity and position estimates using quaternion
            % vkm          = (quat2rotm(quatAttitude)*accel_meas + [0 0 9.81]').*dt + obj.estimate(5:7); 
            vkm          = (quat2rotm(obj.estimate(1:4)')*accel_meas + [0 0 9.81]').*dt + obj.estimate(5:7); 
            rkm          = vkm*dt + obj.estimate(8:10); 
            
            
            quatAttitude = quaternion(obj.estimate(1:4)') * quaternion(quatexp(.5*reshape([0;gyro_meas.*dt],1,4)));
            
            % Covariance Propagation 
            I3 = eye(3); 
            Z3 = zeros(3,3); 
            
            %TODO: check about specific force vs acceleration in F
            F  = [-skewSym(gyro_meas)                          Z3 Z3 -I3 Z3                      ; 
                  -quat2rotm(quatAttitude)*skewSym(accel_meas) Z3 Z3  Z3 -quat2rotm(quatAttitude);
                  Z3                                           I3 Z3  Z3 Z3;
                  Z3                                           Z3 Z3  Z3 Z3;
                  Z3                                           Z3 Z3  Z3 Z3];
            
            STM = eye(15) + F * dt + .5 .* F * F * dt^2;

            % Calculate new covariance
            obj.estimate_covariance = STM * obj.estimate_covariance * STM' + obj.processCovariance(dt);
            
            
            delta = obj.processCovariance(dt); 
            cov = obj.estimate_covariance; 

            obj.estimate(1:4)   = quatAttitude.compact; 
            obj.estimate(5:7)   = vkm; 
            obj.estimate(8:10)  = rkm;             
        end

        function delf = updateWithGravity(obj, acc_meas)
            aPrioriErrorState = zeros(15,1); % Only in attitude

            Z33 = zeros(3,3); 
            I33  = eye(3); 

            H = [skewSym(quat2dcm(obj.estimate(1:4)')*[0;0;-9.81]) , Z33, Z33, Z33, I33];

            delf = (acc_meas - obj.estimate(14:16)) - quat2dcm(obj.estimate(1:4)')*[0;0;-9.81];
            
            Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.accel_cov);

            aPosterioriErrorState = aPrioriErrorState + Kgain*delf;
            
            % aPosterior error state covariance
            obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance; 
            
            obj.estimate(1:4) = quatmultiply(obj.estimate(1:4)', [1 aPosterioriErrorState(1:3)'./2]);

            obj.estimate(11:16) = obj.estimate(11:16) + aPosterioriErrorState(10:15);
        end

        function delm = updateWithMagnetometer(obj, magMeas)
            aPrioriErrorState = zeros(15,1); % Only in attitude

            Z33 = zeros(3,3); 
            I33  = eye(3); 

            NED_mag = [27.5550; -2.4169; -16.0849];
            % NED_mag = [27.5550; 0; -16.0849];

            H = [skewSym(-quat2rotm(obj.estimate(1:4)')*NED_mag), Z33, Z33, Z33, Z33];

            % delm = quat2rotm(obj.estimate(1:4)')'*-NED_mag - magMeas
            delm = magMeas - quat2dcm(obj.estimate(1:4)')*NED_mag;

            Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.mag_cov);

            aPosterioriErrorState = aPrioriErrorState + Kgain*delm;

            % aPosterior error state covariance
            obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance; 

            obj.estimate(1:4) = quatmultiply(obj.estimate(1:4)', [1 aPosterioriErrorState(1:3)'./2]); 

            obj.estimate(11:16) = obj.estimate(11:16) + aPrioriErrorState(10:15); 
        end

        function del_z =  updateWithBarometer(obj, baro_meas)
            aPrioriErrorState = zeros(15,1); % Only in position

            Z13 = zeros(1,3); 

            H = [Z13 Z13 [0 0 1] Z13 Z13]; 

            del_z  = baro_meas - obj.estimate(10);

            Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.barometer_variance);

            aPosterioriErrorState = aPrioriErrorState + Kgain*del_z; 
            % aPosterioriErrorState(9)

            obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance; 

            obj.estimate(5:end) = obj.estimate(5:end) + aPosterioriErrorState(4:end); 

        end
        
        function del_v = updateWithGPS(obj, vel_measurement, pos_measurement)
            assert(~any(size(pos_measurement)  ~= [3,1]))
            assert(~any(size(vel_measurement) ~= [3,1]))
            
            aPrioriErrorState = zeros(15,1); % Only in attitude

            Z33 = zeros(3,3); 
            I3  = eye(3);

            H = [Z33 I3 Z33 Z33 Z33;
                 Z33 Z33 I3 Z33 Z33]; 

            del_vx = [vel_measurement; pos_measurement] - obj.estimate(5:10);

            Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.gpsVariance);

            aPosterioriErrorState = aPrioriErrorState + Kgain*del_vx; 

            obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance; 

            % Make updates only to position and velocity
            obj.estimate(5:10) = obj.estimate(5:10) + aPosterioriErrorState(4:9); 
            
            obj.estimate(14:16) = obj.estimate(14:16) + aPosterioriErrorState(13:15);
        end

        function [state, covariance] = getFilterState(obj)
            state       = obj.estimate; 
            covariance  = diag(obj.estimate_covariance); 
        end

    end
end
