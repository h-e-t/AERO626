classdef EGMF_CLASS < handle
    % Extended Gaussian Mixture Filter (EGMF)
    %
    % This class maintains a mixture of MEKF_CLASS filters.
    % Each mixture component has: 
    %   - a weight w_i
    %   - a state estimate (inside MEKF)
    %   - a covariance estimate (inside MEKF)
    %
    % All propagate() and update() calls simply loop across
    % each MEKF component and update its weight using the
    % measurement likelihood.
    properties
        components           % Array of MEKF_CLASS objects
        weights              % Nx1 mixture weights
        numComponents  

        % pruning & merging parameters
        minWeight = 1e-4
        maxComponents = 10
    end
    
    methods(Access=public)
        %% ---------------------------------------------------------------
        function obj = EGMF_CLASS(initGuess)
            arguments
                initGuess MEKF_CLASS
            end

            N = 2 * 15 + 1;

            obj.components       = cell(1,N); 

            errorStateCovariance = diag(initGuess.estimate_covariance); 

            for idx = 1:N
                obj.components{idx} = initGuess.clone(); 
            end

            % Component 1 remains the original component
            obj.components{2}.rotate([errorStateCovariance(1) 0 0]);
            obj.components{3}.rotate([0 errorStateCovariance(2) 0]);
            obj.components{4}.rotate([ 0 0 errorStateCovariance(3)]); 

            obj.components{5}.rotate([-errorStateCovariance(1) 0 0]); 
            obj.components{6}.rotate([0 -errorStateCovariance(2) 0]); 
            obj.components{7}.rotate([ 0 0 -errorStateCovariance(3)]);

            for idx = 8:2:31
                covIdx = idx/2; 

                obj.components{idx}.estimate(covIdx+1)   =   obj.components{idx}.estimate(covIdx+1) + errorStateCovariance(covIdx); 
                obj.components{idx+1}.estimate(covIdx+1) = obj.components{idx+1}.estimate(covIdx+1) - errorStateCovariance(covIdx); 
            end

            obj.weights             = ones(1,N)/N;
            obj.numComponents       = N;
        end

        %% ---------------------------------------------------------------
        function propagate(obj, gyro_meas, accel_meas, dt)
            % Propagate each mixture component independently
            for i = 1:obj.numComponents
                obj.components{i}.propagate(gyro_meas, accel_meas, dt);
            end
        end

        function updateWithGravity(obj, acc_meas)
            likelihoods = zeros(obj.numComponents,1);

            for i = 1:obj.numComponents
                mekf = obj.components{i};

                % Predict gravity in acc frame
                R = quat2rotm(mekf.estimate(1:4)');
                g_est = R * [0;0;-9.81];

                % Likelihood of measurement
                innov = acc_meas - g_est;
                S = mekf.accel_cov * eye(3);
                likelihoods(i) = obj.gaussianLikelihood(innov, S);

                % Execute MEKF gravity update
                mekf.updateWithGravity(acc_meas);
            end

            obj.updateWeights(likelihoods);
        end

        %% ---------------------------------------------------------------
        function updateWithMagnetometer(obj, mag_meas)

            likelihoods = zeros(obj.numComponents,1);

            for i = 1:obj.numComponents
                mekf = obj.components{i};

                % Compute predicted mag body measurement
                pred_mag = quat2rotm(mekf.estimate(1:4)') * mag_meas;

                innov = mag_meas - pred_mag;
                S = mekf.mag_cov * eye(3);

                likelihoods(i) = obj.gaussianLikelihood(innov, S);

                % MEKF internal update
                mekf.updateWithMagnetometer(mag_meas);
            end

            obj.updateWeights(likelihoods);
        end
        
        %% ---------------------------------------------------------------
        function updateWithBarometer(obj, baro_meas)

            likelihoods = zeros(obj.numComponents,1);

            for i = 1:obj.numComponents
                mekf = obj.components{i};

                z_pred = mekf.estimate(10);  % predicted altitude
                innov  = baro_meas - z_pred;
                S = mekf.barometer_variance; 

                likelihoods(i) = obj.gaussianLikelihood(innov, S);

                mekf.updateWithBarometer(baro_meas);
            end

            obj.updateWeights(likelihoods);
        end

        %% ---------------------------------------------------------------
        function updateWithGPS(obj, vel_meas, pos_meas)

            likelihoods = zeros(obj.numComponents,1);

            z = [vel_meas; pos_meas];

            for i = 1:obj.numComponents
                mekf = obj.components{i};

                z_pred = mekf.estimate(5:10); % predicted vel + pos
                innov  = z - z_pred;

                S = mekf.gpsVariance;  % 6×6 covariance

                likelihoods(i) = obj.gaussianLikelihood(innov, S);

                mekf.updateWithGPS(vel_meas, pos_meas);
            end

            obj.updateWeights(likelihoods);
        end

        %% ---------------------------------------------------------------
        function updateWeights(obj, likelihoods)
            % Bayesian weight update
            obj.weights = obj.weights .* likelihoods;
            w_sum = sum(obj.weights);

            if w_sum < 1e-12
                % Degeneracy — reset to uniform
                obj.weights = ones(obj.numComponents,1)/obj.numComponents;
            else
                obj.weights = obj.weights / w_sum;
            end

            obj.pruneSmallWeights();
        end

        %% ---------------------------------------------------------------
        function pruneSmallWeights(obj)
            keep = obj.weights > obj.minWeight;
            obj.components = obj.components(keep);
            obj.weights = obj.weights(keep);
            obj.numComponents = length(obj.weights);

            % Normalize again after pruning
            obj.weights = obj.weights/sum(obj.weights);
        end

        %% ---------------------------------------------------------------
        function [state, cov] = getFilterState(obj)
            % Weighted mean state
            state = zeros(16,1);
            for i = 1:obj.numComponents
                state = state + obj.weights(i) * obj.components{i}.estimate;
            end

            % Weighted covariance approximation
            cov = zeros(15,15);
            for i = 1:obj.numComponents
                diff = obj.components{i}.estimate(1:15) - state(1:15);
                cov = cov + obj.weights(i) * ...
                      (obj.components{i}.estimate_covariance + diff*diff');
            end

            cov = diag(cov); 
        end

    end

    methods(Access=private)
        function L = gaussianLikelihood(obj, innovation, S)
            k = length(innovation);
            L = exp(-0.5 * innovation' / S * innovation) / sqrt((2*pi)^k * det(S));
        end
    end

end