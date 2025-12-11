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

        mag_cov
        gpsVariance

        consecutive_grav_updates
        covarianceFloor
    end

    methods
        
        %==================================================================
        % State Vector Definitions
        %
        % TRUE STATE VECTOR:
        %  [1:4]   quaternion attitude
        %  [5:7]   velocity 
        %  [8:10]  position
        % [11:13]  gyro bias
        % [14:16]  accel bias
        %
        % ERROR STATE VECTOR:
        %  [1:3]   attitude error (3-vector)
        %  [4:6]   velocity error
        %  [7:9]   position error
        % [10:12]  gyro bias error
        % [13:15]  accel bias error
        %==================================================================

        function obj = MEKF_CLASS(init_estimate, estimate_covariance, ...
                                  gyro_cov, gyro_bias_cov, ...
                                  accel_cov, accel_bias_cov, ...
                                  mag_cov, gps_cov, baro_cov)

            obj.estimate            = init_estimate;
            obj.estimate_covariance = estimate_covariance;

            obj.gyro_cov            = gyro_cov;
            obj.gyro_bias           = zeros(3, 1);
            obj.gyro_bias_cov       = gyro_bias_cov;

            obj.accel_cov           = accel_cov;
            obj.accel_bias          = zeros(3, 1);
            obj.accel_bias_cov      = accel_bias_cov;

            obj.mag_cov             = mag_cov;

            obj.barometer_variance  = baro_cov;
            obj.gpsVariance         = gps_cov;

            obj.consecutive_grav_updates = 0;

            % Floor applied to covariance diagonal entries
            obj.covarianceFloor = [ ...
                0.0012 0.0012 0.0012 ...  % attitude
                0.05   0.05   0.05   ...  % velocity
                0.2    0.2    1.0    ...  % position 
                0.05   0.05   0.05   ...  % gyro bias
                0.05   0.05   0.05 ]';    % accel bias
        end

        %==================================================================
        % Process Covariance
        %==================================================================
        function Qd = processCovariance(obj, dt)
            Z33 = zeros(3,3);

            sw  = obj.gyro_cov^2;
            sf  = obj.accel_cov^2;

            sbw = obj.gyro_bias_cov^2;
            sbf = obj.accel_bias_cov^2;

            Qd = [ ...
                (sbw*dt^3)/3 + sw*dt,             Z33,                        Z33, -(dt^2*sbw)/2,    Z33;
                Z33,           (sbf*dt^3)/3 + sf*dt, (sbf*dt^4)/8 + (sf*dt^2)/2,  Z33, -(dt^2*sbf)/2;
                Z33, (sbf*dt^4)/8 + (sf*dt^2)/2, (sbf*dt^5)/20 + (sf*dt^3)/3,    Z33, -(dt^3*sbf)/6;
                -(dt^2*sbw)/2, Z33,               Z33,                         dt*sbw,              Z33;
                Z33,           -(dt^2*sbf)/2,    -(dt^3*sbf)/6,                Z33, dt*sbf ];
        end

        %==================================================================
        % State Propagation
        %==================================================================
        function propagate(obj, gyro_meas, acc_meas, dt)
            
            gyro_meas  = gyro_meas - obj.gyro_bias;
            accel_meas = acc_meas  - obj.accel_bias;

            assert(any(size(gyro_meas) == [3,1]))
            assert(any(size(accel_meas) == [3,1]))

            % Quaternion attitude update
            quatAttitude = quaternion(obj.estimate(1:4)') * ...
                           quaternion(quatexp(.5 * reshape([0; gyro_meas.*dt], 1, 4)));

            % Rotation matrix based on prior estimate
            R = quat2rotm(obj.estimate(1:4)');
            g = [0; 0; 9.81];

            % Velocity and position propagation
            vkm = (R * accel_meas + g) * dt + obj.estimate(5:7);
            rkm = vkm * dt + obj.estimate(8:10);

            I3 = eye(3);
            Z3 = zeros(3,3);

            F = [ ...
                -skewSym(gyro_meas),                          Z3, Z3, -I3,              Z3;
                -quat2rotm(quatAttitude) * skewSym(accel_meas), Z3, Z3,  Z3, -quat2rotm(quatAttitude);
                 Z3,                                          I3, Z3,  Z3,              Z3;
                 Z3,                                          Z3, Z3,  Z3,              Z3;
                 Z3,                                          Z3, Z3,  Z3,              Z3 ];

            STM = eye(15) + F*dt + 0.5*(F*F)*dt^2;

            % Covariance propagation
            obj.estimate_covariance = STM * obj.estimate_covariance * STM' + ...
                                      obj.processCovariance(dt);

            % Write state back
            obj.estimate(1:4)   = quatAttitude.compact;
            obj.estimate(5:7)   = vkm;
            obj.estimate(8:10)  = rkm;

            % Normalize quaternion
            obj.estimate(1:4) = obj.estimate(1:4) ./ norm(obj.estimate(1:4));
        end

        %==================================================================
        % Gravity-Based Update
        %==================================================================
        function updateWithGravity(obj, acc_meas)

            N = 15;
            isValidReading = abs(9.81 - norm(acc_meas)) < 0.4;

            if isValidReading
                obj.consecutive_grav_updates = obj.consecutive_grav_updates + 1;
            else
                obj.consecutive_grav_updates = 0;
                return
            end

            if obj.consecutive_grav_updates < N
                return
            end

            aPrioriErrorState = zeros(15,1);

            Z33 = zeros(3,3);
            I33 = eye(3);

            H = [ skewSym(quat2dcm(obj.estimate(1:4)') * [0;0;-9.81]), ...
                  Z33, Z33, Z33, I33 ];

            delf = (acc_meas - obj.estimate(14:16)) - ...
                   quat2dcm(obj.estimate(1:4)') * [0;0;-9.81];

            Kgain = obj.estimate_covariance * H' / ...
                    (H * obj.estimate_covariance * H' + obj.accel_cov);

            aPosterioriErrorState = aPrioriErrorState + Kgain * delf;

            obj.estimate_covariance = (eye(15) - Kgain*H) * obj.estimate_covariance;

            obj.estimate(1:4) = quatmultiply(obj.estimate(1:4)', ...
                                             [1 0 aPosterioriErrorState(2:3)'/2]);

            obj.estimate(11:16) = obj.estimate(11:16) + aPosterioriErrorState(10:15);
        end

        %==================================================================
        % Magnetometer Update
        %==================================================================
        function updateWithMagnetometer(obj, magMeas)

            NED_mag         = [27.5550; -2.4169; -16.0849];
            NED_declination = -0.0875;

            Estimated_Mag = quat2rotm(obj.estimate(1:4)') * magMeas;
            EstimatedDeclination = atan2(Estimated_Mag(2), Estimated_Mag(1));

            del_beta = [0; 0; (NED_declination - EstimatedDeclination)];

            H = [0 0 0;
                 0 0 0;
                 0 0 1];

            Kgain = obj.estimate_covariance(1:3,1:3) * H / ...
                    (obj.mag_cov + H * obj.estimate_covariance(1:3,1:3) * H');

            obj.estimate(1:4) = quatmultiply(obj.estimate(1:4)', ...
                                             [1 (Kgain*del_beta)'/2]);

            obj.estimate_covariance(1:3,1:3) = (eye(3) - Kgain) * ...
                                               obj.estimate_covariance(1:3,1:3) * 0.98;
        end

        %==================================================================
        % Barometer Update
        %==================================================================
        function updateWithBarometer(obj, baro_meas)

            aPrioriErrorState = zeros(15,1);
            Z13 = zeros(1,3);

            H = [Z13 Z13 [0 0 1] Z13 Z13];

            del_z = baro_meas - obj.estimate(10);

            Kgain = obj.estimate_covariance * H' / ...
                    (H * obj.estimate_covariance * H' + obj.barometer_variance);

            aPosterioriErrorState = aPrioriErrorState + Kgain * del_z;

            obj.estimate_covariance = (eye(15) - Kgain*H) * obj.estimate_covariance;

            obj.estimate(5:end) = obj.estimate(5:end) + aPosterioriErrorState(4:end);

            obj.checkFloor();
        end

        %==================================================================
        % GPS Update
        %==================================================================
        function updateWithGPS(obj, vel_measurement, pos_measurement)

            assert(~any(size(pos_measurement)  ~= [3,1]))
            assert(~any(size(vel_measurement) ~= [3,1]))

            aPrioriErrorState = zeros(15,1);

            Z33 = zeros(3,3);
            I3  = eye(3);

            H = [Z33 I3 Z33 Z33 Z33;
                 Z33 Z33 I3 Z33 Z33];

            del_vx = [vel_measurement; pos_measurement] - obj.estimate(5:10);

            S = H * obj.estimate_covariance * H' + obj.gpsVariance;
            K = obj.estimate_covariance * H' / S;

            aPosterioriErrorState = aPrioriErrorState + K * del_vx;

            I = eye(15);
            obj.estimate_covariance = (I - K*H) * obj.estimate_covariance * (I - K*H)' + ...
                                      K * obj.gpsVariance * K';

            obj.estimate(5:10) = obj.estimate(5:10) + aPosterioriErrorState(4:9);

            obj.checkFloor();
        end

        %==================================================================
        % Utility Functions
        %==================================================================
        function checkFloor(obj)
            d = diag(obj.estimate_covariance);
            d = max(d, obj.covarianceFloor);
            obj.estimate_covariance(1:16:end) = d;
        end

        function [state, covariance] = getFilterState(obj)
            state      = obj.estimate;
            covariance = diag(obj.estimate_covariance);
        end

        function newObj = clone(obj)
            newObj = MEKF_CLASS(obj.estimate, obj.estimate_covariance, ...
                                obj.gyro_cov, obj.gyro_bias_cov, ...
                                obj.accel_cov, obj.accel_bias_cov, ...
                                obj.mag_cov, obj.gpsVariance, obj.barometer_variance);

            newObj.gyro_bias  = obj.gyro_bias;
            newObj.accel_bias = obj.accel_bias;
            newObj.consecutive_grav_updates = obj.consecutive_grav_updates;
            newObj.covarianceFloor = obj.covarianceFloor;
        end

        function rotate(obj, rotation)
            arguments
                obj
                rotation
            end

            assert(any(size(rotation) ~= [3 1]))

            quatAttitude = quaternion(obj.estimate(1:4)') * quaternion(quatexp(.5*[0 rotation]));

            obj.estimate(1:4)   = quatAttitude.compact; 
        end
    end
end


% classdef MEKF_CLASS < handle
% 
%     properties
%         estimate
%         estimate_covariance
% 
%         gyro_cov
%         gyro_bias
%         gyro_bias_cov
% 
%         accel_cov
%         accel_bias
%         accel_bias_cov
% 
%         barometer_variance
%         barometer_bias
% 
%         mag_cov
% 
%         gpsVariance
% 
%         consecutive_grav_updates
% 
%         covarianceFloor
%     end
% 
%     methods
%         % STATE_VECTOR 
%         %  [1:4]  quaternion attitude
%         %  [5:7]  velocity 
%         %  [8:10] position
%         % [11:13] gyro  bias
%         % [14:16] accel bias
% 
%         % ERROR STATE_VECTOR 
%         %  [1:3]  quaternion attitude
%         %  [4:6]  velocity 
%         %  [7:9]  position
%         % [10:12] gyro  bias
%         % [13:15] accel bias
%         function obj = MEKF_CLASS(init_estimate, estimate_covariance, ...
%                                   gyro_cov,      gyro_bias_cov,  ...
%                                   accel_cov,     accel_bias_cov, ...
%                                   mag_cov, ...
%                                   gps_cov, ...
%                                   baro_cov)
% 
%             obj.estimate            = init_estimate;
%             obj.estimate_covariance = estimate_covariance ; 
% 
%             obj.gyro_cov            = gyro_cov;
%             obj.gyro_bias           = zeros(3, 1);
%             obj.gyro_bias_cov       = gyro_bias_cov; 
% 
%             obj.accel_cov           = accel_cov; 
%             obj.accel_bias          = zeros(3,1); 
%             obj.accel_bias_cov      = accel_bias_cov; 
% 
%             obj.mag_cov             = mag_cov; 
% 
%             obj.barometer_variance  = baro_cov; 
% 
%             obj.gpsVariance         = gps_cov; 
% 
%             obj.consecutive_grav_updates  = 0;
% 
%             obj.covarianceFloor = [0.0012 0.0012 0.0012 .05 .05 .05 .2 .2 1 .05 .05 .05 .05 .05 .05]';
%         end
% 
%         function Qd = processCovariance(obj, dt)
%             Z33 = zeros(3,3); 
% 
%             sw = obj.gyro_cov^2; 
%             sf = obj.accel_cov^2; 
% 
%             sbw = obj.gyro_bias_cov^2; 
%             sbf = obj.accel_bias_cov ^ 2; 
% 
%             Qd = ... 
%             [(sbw*dt^3)/3 + sw*dt,                         Z33,                        Z33, -(dt^2*sbw)/2,           Z33;
%                               Z33,       (sbf*dt^3)/3 + sf*dt,  (sbf*dt^4)/8 + (sf*dt^2)/2,           Z33, -(dt^2*sbf)/2;
%                               Z33, (sbf*dt^4)/8 + (sf*dt^2)/2, (sbf*dt^5)/20 + (sf*dt^3)/3 ,           Z33, -(dt^3*sbf)/6;
%                    -(dt^2*sbw)/2,                         Z33,                         Z33,        dt*sbw,           Z33;
%                               Z33,              -(dt^2*sbf)/2,               -(dt^3*sbf)/6,           Z33,        dt*sbf];
%         end
% 
%         function propagate(obj, gyro_meas, acc_meas, dt)
%             gyro_meas   = gyro_meas - obj.gyro_bias;
%             accel_meas  = acc_meas  - obj.accel_bias;
% 
%             assert(any(size(gyro_meas) ==[3,1]))
%             assert(any(size(accel_meas)==[3,1]))
% 
% 
% 
%             % Propagate quaternion attiude using measurement
%             quatAttitude = quaternion(obj.estimate(1:4)') * quaternion(quatexp(.5*reshape([0;gyro_meas.*dt],1,4)));
% 
%             R = quat2rotm(obj.estimate(1:4)'); %
%             g = [0 0 9.81]'; 
% 
%             % update the velocity and position estimates using quaternion
%             vkm          = (R*accel_meas + g)*dt + obj.estimate(5:7); 
%             rkm          = vkm*dt + obj.estimate(8:10); 
% 
% 
%             I3 = eye(3); 
%             Z3 = zeros(3,3); 
% 
%             F  = [-skewSym(gyro_meas)                          Z3 Z3 -I3 Z3                      ; 
%                   -quat2rotm(quatAttitude)*skewSym(accel_meas) Z3 Z3  Z3 -quat2rotm(quatAttitude);
%                   Z3                                           I3 Z3  Z3 Z3;
%                   Z3                                           Z3 Z3  Z3 Z3;
%                   Z3                                           Z3 Z3  Z3 Z3];
% 
%             STM = eye(15) + F * dt + .5 .* F * F * dt^2;
% 
%             % Calculate new covariance
%             obj.estimate_covariance = STM * obj.estimate_covariance * STM' + obj.processCovariance(dt);
% 
%             obj.estimate(1:4)   = quatAttitude.compact; 
%             obj.estimate(5:7)   = vkm; 
%             obj.estimate(8:10)  = rkm;    
% 
%             if norm(obj.estimate(1:4)') ~= 1
%                 obj.estimate(1:4) = obj.estimate(1:4) / norm(obj.estimate(1:4)');
%             end
%         end
% 
%         function updateWithGravity(obj, acc_meas)
%             % Basic Outlier Rejection 
%             N = 15; 
%             isValidReading = abs(9.81 - norm(acc_meas)) < .4; 
% 
%             if  isValidReading 
%                 obj.consecutive_grav_updates = obj.consecutive_grav_updates + 1; 
%             else
%                 obj.consecutive_grav_updates = 0;
%                 return
%             end
% 
%             if obj.consecutive_grav_updates < N
%                 return
%             end
% 
%             aPrioriErrorState = zeros(15,1); % Only in attitude
% 
%             Z33 = zeros(3,3); 
%             I33  = eye(3); 
% 
%             H = [skewSym(quat2dcm(obj.estimate(1:4)')*[0;0;-9.81]) , Z33, Z33, Z33, I33];
% 
%             delf = (acc_meas - obj.estimate(14:16)) - quat2dcm(obj.estimate(1:4)')*[0;0;-9.81];
% 
%             Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.accel_cov);
% 
%             aPosterioriErrorState = aPrioriErrorState + Kgain*delf;
% 
%             % aPosterior error state covariance
%             obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance; 
% 
%             obj.estimate(1:4) = quatmultiply(obj.estimate(1:4)', [1 0 aPosterioriErrorState(2:3)'./2]);
% 
%             obj.estimate(11:16) = obj.estimate(11:16) + aPosterioriErrorState(10:15);
%         end
% 
%         function updateWithMagnetometer(obj, magMeas)
%             NED_mag         = [27.5550; -2.4169; -16.0849];
% 
%             NED_declination = -0.0875; 
% 
%             Estimated_Mag           = quat2rotm(obj.estimate(1:4)')*magMeas;
% 
%             EstimatedDeclination    = atan2(Estimated_Mag(2), Estimated_Mag(1));
% 
%             del_beta = [0 0 (NED_declination - EstimatedDeclination)]';
% 
%             H = [0 0 0;
%                  0 0 0;
%                  0 0 1];
% 
%             Kgain    = obj.estimate_covariance(1:3,1:3) * H/(obj.mag_cov + H * obj.estimate_covariance(1:3,1:3) * H');
% 
%             obj.estimate(1:4) = quatmultiply(obj.estimate(1:4)', [1 (Kgain*del_beta)'./2]); 
% 
%             obj.estimate_covariance(1:3, 1:3) = (eye(3) - Kgain)*obj.estimate_covariance(1:3, 1:3)*.98; 
%         end
% 
%         function updateWithBarometer(obj, baro_meas)
%             aPrioriErrorState = zeros(15,1); 
% 
%             Z13 = zeros(1,3); 
% 
%             H = [Z13 Z13 [0 0 1] Z13 Z13]; 
% 
%             del_z  = baro_meas - obj.estimate(10);
% 
%             Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.barometer_variance);
% 
%             aPosterioriErrorState = aPrioriErrorState + Kgain*del_z; 
% 
%             obj.estimate_covariance = (eye(15) - Kgain*H)*obj.estimate_covariance;
% 
%             obj.estimate(5:end) = obj.estimate(5:end) + aPosterioriErrorState(4:end); 
% 
%             obj.checkFloor()
%         end
% 
%         function updateWithGPS(obj, vel_measurement, pos_measurement)
%             assert(~any(size(pos_measurement)  ~= [3,1]))
%             assert(~any(size(vel_measurement) ~= [3,1]))
% 
%             aPrioriErrorState = zeros(15,1);
% 
%             Z33 = zeros(3,3); 
%             I3  = eye(3);
% 
%             H = [Z33 I3 Z33 Z33 Z33;
%                  Z33 Z33 I3 Z33 Z33]; 
% 
%             del_vx = [vel_measurement; pos_measurement] - obj.estimate(5:10);
% 
%             Kgain = obj.estimate_covariance*H'/(H * obj.estimate_covariance * H' + obj.gpsVariance);
% 
%             aPosterioriErrorState = aPrioriErrorState + Kgain*del_vx;
% 
% 
%             I = eye(15);
%             S = H * obj.estimate_covariance * H' + obj.gpsVariance;
%             K = obj.estimate_covariance * H' / S;
% 
%             obj.estimate_covariance = (I - K*H)*obj.estimate_covariance*(I - K*H)' + K*obj.gpsVariance*K'; 
%             obj.estimate(5:10) = obj.estimate(5:10) + aPosterioriErrorState(4:9);
% 
% 
%             obj.checkFloor()
%         end
% 
%         function checkFloor(obj)
%             d = diag(obj.estimate_covariance); 
%             d = max(d, obj.covarianceFloor); 
%             obj.estimate_covariance(1:16:end) = d; 
%         end
% 
%         function [state, covariance] = getFilterState(obj)
%             state       = obj.estimate; 
%             covariance  = diag(obj.estimate_covariance); 
%         end
% 
%         % EGMF Helper Functions 
% 
%         function newObj = clone(obj)
%             % Create a new MEKF with identical internal state
%             newObj = MEKF_CLASS( ...
%                 obj.estimate, ...
%                 obj.estimate_covariance, ...
%                 obj.gyro_cov, ...
%                 obj.gyro_bias_cov, ...
%                 obj.accel_cov, ...
%                 obj.accel_bias_cov, ...
%                 obj.mag_cov, ...
%                 obj.gpsVariance, ...
%                 obj.barometer_variance );
% 
%             % Copy all dynamic internal variables
%             newObj.gyro_bias  = obj.gyro_bias;
%             newObj.accel_bias = obj.accel_bias;
%             newObj.consecutive_grav_updates = obj.consecutive_grav_updates;
%             newObj.covarianceFloor = obj.covarianceFloor;
%         end
% 
%         function rotate(obj, rotation)
%             arguments
%                 obj
%                 rotation 
%             end
% 
%             assert(any(size(rotation) ~= [3 1]))
% 
%             quatAttitude = quaternion(obj.estimate(1:4)') * quaternion(quatexp(.5*[0 rotation]));
% 
%             obj.estimate(1:4)   = quatAttitude.compact; 
%         end
%     end
% end
% 
% 
