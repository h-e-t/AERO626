classdef ExtendedGaussianMixtureFilter < handle
    % This class implements the Gaussian mixture filter in the linear, 
    % extended, and unscented form for an arbitrary state size and 
    % arbitrary state dynamics

    properties
        means       (:,1,:)
        covariances (:,:,:)
        weights     (1,:)
        stateSize   (1,1)       
        numWeights  (1,1)

        alpha      
        beta   
        kappa  

        npl
        lambda

        unscentedMeanWeights
        unscentedCovarianceWeights
    end

    methods
        function obj = ExtendedGaussianMixtureFilter(means,covariances,weights, alpha, beta, kappa)            
            arguments
                means 
                covariances 
                weights 
                alpha   = 0
                beta    = 0
                kappa   = 0
            end
            
            obj.means           = means; 
            obj.covariances     = covariances; 
            obj.weights         = weights; 
            obj.stateSize       = length(means(:,1)); 
            obj.numWeights      = length(weights); 

            obj.alpha           = alpha;
            obj.beta            = beta;
            obj.kappa           = kappa;
            
            obj.lambda          = alpha^2*(obj.stateSize + kappa) - obj.stateSize; 
            obj.npl             = obj.stateSize + obj.lambda; 
            
            numSigmaPoints      = 1 + 2*obj.stateSize;
            
            obj.unscentedMeanWeights            = obj.lambda / obj.npl * ones(1,numSigmaPoints);
            obj.unscentedMeanWeights(1,2:end)   = obj.unscentedMeanWeights(1,2:end)/(2*obj.lambda);
            
            obj.unscentedCovarianceWeights      = obj.unscentedMeanWeights; 
            obj.unscentedCovarianceWeights(1)   = obj.unscentedCovarianceWeights(1) + (1 - obj.alpha^2 + obj.beta);
        end

        function predictPDF_Linear(obj, F, Pww)
            % For a given linear state transition model and process noise,
            % this function updates the state of the filter with a
            % propagation. 

            newMeans = zeros(size(obj.means)); 
            newCovariances = zeros(size(obj.covariances));

            for idx = 1:length(obj.weights)
                newMeans(:,:,idx)       = F * obj.means(:,:,idx); 
                newCovariances(:,:,idx) = F * obj.covariances(:,:,idx) * F' + Pww;  
            end
            
            obj.means       = newMeans; 
            obj.covariances = newCovariances; 
        end 

        function predictPDF_Extended(obj, f, Fx, Pww)
            % For a given non linear state transition function and process 
            % noise. This function updates the state of the filter with a
            % propagation. 

            newMeans = zeros(size(obj.means)); 
            newCovariances = zeros(size(obj.covariances));

            for idx = 1:obj.numWeights
                newMeans(:,:,idx)       = f(obj.means(:,:,idx)); 
                newCovariances(:,:,idx) = Fx(obj.means(:,:,idx)) * obj.covariances(:,:,idx) * Fx(obj.means(:,:,idx))' + Pww;  
            end
            
            obj.means       = newMeans; 
            obj.covariances = newCovariances; 
        end

        function predictPDF_Unscented(obj, f, Pww)
            numSigmaPoints      = 1 + 2*obj.stateSize;

            for idx = 1:obj.numWeights
                sigmaPoints = zeros(obj.stateSize, numSigmaPoints); 
                
                % Generate Sigma Points
                SxxL = chol(obj.covariances(:,:,idx), "lower");
                                
                sigmaPoints(:,1) = obj.means(:,1,idx);
                
                sigmaPointIdx = 2; 
                for jdx = 1:obj.stateSize
                    sigmaPoints(:,sigmaPointIdx)    = obj.means(:,1,idx) + sqrt(obj.npl)*SxxL(:,jdx); 
                    sigmaPoints(:,sigmaPointIdx+1)  = obj.means(:,1,idx) - sqrt(obj.npl)*SxxL(:,jdx);
                    
                    sigmaPointIdx = sigmaPointIdx+1; 
                end
                
                % Propagate Sigma Points
                for jdx = 1:numSigmaPoints
                    sigmaPoints(:,jdx) = f(sigmaPoints(:,jdx));
                end 
                
                % Propagate Mean 
                obj.means(:,1, idx) = sum(sigmaPoints.*obj.unscentedMeanWeights);

                % Propagate Covariance 
                obj.covariances(:,:,idx) = Pww;
                for jdx = 1:numSigmaPoints
                    obj.covariances(:,:,idx) = obj.covariances(:,:,idx) + ...
                        obj.unscentedCovarianceWeights(jdx).*(sigmaPoints(:,jdx) - obj.means(:,1, idx))...
                        * (sigmaPoints(:,jdx) - obj.means(:,1, idx))';
                end

            end
            
        end
        
        function updatePDF_Linear(obj, H, Pvv, Z)
            % For a given measurement, measurement model, and measurement
            % uncertainty, this function updates the internal state of the
            % system

            littlek                     = zeros(1,obj.numWeights);
            aPosteriorWeights           = zeros(1,obj.numWeights); 
            aPosteriorCovariances       = zeros(size(obj.covariances)); 
            aPosteriorMeans             = zeros(size(obj.means)); 
            bigK                        = zeros(obj.stateSize, length(Z), obj.numWeights); 


            for idx = 1:obj.numWeights
                littlek(idx)                    = mvnpdf(Z, H*obj.means(:,idx), H * obj.covariances(:,:,idx) * H' + Pvv); 
                bigK(:,:,idx)                   = obj.covariances(:,:,idx) * H' / (H * obj.covariances(:,:,idx) * H' + Pvv);

                aPosteriorMeans(:,:,idx)        = obj.means(:,:,idx) + bigK(:,:,idx)*(Z - H * obj.means(:,:,idx)); 
                aPosteriorCovariances(:,:,idx)  = obj.covariances(:,:,idx) - bigK(:,:,idx) * H * obj.covariances(:,:,idx); 
            end
            
            aPosteriorWeights(1,:)  = littlek.*obj.weights/(sum(littlek.*obj.weights)); 
            
            obj.means               = aPosteriorMeans; 
            obj.covariances         = aPosteriorCovariances; 
            obj.weights             = aPosteriorWeights;  
            
        end

        function updatePDF_Extended(obj,h, Hx, Pvv, Z)
            % For a given measurement, measurement model, and measurement
            % uncertainty, this function updates the internal state of the
            % system
            
            littlek                     = zeros(1,obj.numWeights);
            aPosteriorWeights           = zeros(1,obj.numWeights); 
            aPosteriorCovariances       = zeros(size(obj.covariances)); 
            aPosteriorMeans             = zeros(size(obj.means)); 
            bigK                        = zeros(obj.stateSize, length(Z), obj.numWeights); 
            
            
            for idx = 1:obj.numWeights
                Pxz = obj.covariances(:,:,idx) * Hx(obj.means(:,idx))'; 
                Pzz = Hx(obj.means(:,idx)) * obj.covariances(:,:,idx) * Hx(obj.means(:,idx))' + Pvv;
                
                % if obj.numWeights ~= 1
                littlek(idx)                = mvnpdf(Z, h(obj.means(:,idx)), Pzz); 
                % else
                %     littlek(idx)                = normpdf(Z, h(obj.means(:,idx)), Pzz);
                % end

                bigK(:,:,idx)                   = Pxz / Pzz;
                
                aPosteriorMeans(:,:,idx)        = obj.means(:,:,idx) + bigK(:,:,idx)*(Z - h(obj.means(:,:,idx))); 
                aPosteriorCovariances(:,:,idx)  = obj.covariances(:,:,idx)... 
                                                    - Pxz*bigK(:,:,idx)' ...
                                                    - bigK(:,:,idx) * Pxz' ...
                                                    + bigK(:,:,idx) * Pzz * bigK(:,:,idx)'; 
            end
            
            aPosteriorWeights(1,:)  = littlek.*obj.weights/(sum(littlek.*obj.weights)); 
            
            obj.means               = aPosteriorMeans; 
            obj.covariances         = aPosteriorCovariances; 
            obj.weights             = aPosteriorWeights;  

        end

        function updatePDF_Unscented(obj, h, Pvv, Z)
            numSigmaPoints      = 1 + 2*obj.stateSize;

            for idx = 1:obj.numWeights
                sigmaPoints = zeros(obj.stateSize, numSigmaPoints); 
                measurementSigmaPoints = zeros(obj.stateSize, numSigmaPoints);

                % Generate Sigma Points
                SxxL = chol(obj.covariances(:,:,idx), "lower");
                                
                sigmaPoints(:,1) = obj.means(:,1,idx);
                
                sigmaPointIdx = 2; 
                for jdx = 1:obj.stateSize
                    sigmaPoints(:,sigmaPointIdx)    = obj.means(:,1,idx) + sqrt(obj.npl)*SxxL(:,jdx); 
                    sigmaPoints(:,sigmaPointIdx+1)  = obj.means(:,1,idx) - sqrt(obj.npl)*SxxL(:,jdx);
                    
                    sigmaPointIdx = sigmaPointIdx+1; 
                end
                
                % Transform Sigma Points to corresponding measurement
                for jdx = 1:numSigmaPoints
                    measurementSigmaPoints(:,jdx) = h(sigmaPoints(:,jdx));
                end 
                
                % Obtain Transformed Mean (Measurement Transformation)
                mzkm = sum(measurementSigmaPoints.*obj.unscentedMeanWeights);
                
                % Obtain measurement covariance
                Pzzk = Pvv;
                for jdx = 1:numSigmaPoints
                    Pzzk = Pzzk + ...
                        obj.unscentedCovarianceWeights(jdx)...
                        .*(measurementSigmaPoints(:,jdx) - mzkm)...
                        * (measurementSigmaPoints(:,jdx) - mzkm)';
                end

                % Obtain cross covariance
                Pxzk = zeros(size(Pvv));
                for jdx = 1:numSigmaPoints
                    Pxzk = Pxzk + ...
                        obj.unscentedCovarianceWeights(jdx)...
                        .*(sigmaPoints(:,jdx) - obj.means(:,1,idx))...
                        * (measurementSigmaPoints(:,jdx) - mzkm)';
                end

                Kgain = Pxzk/Pzzk; 

                obj.means(:,1,idx)       = obj.means(:,1,idx)  ...
                                           + Kgain*(Z - mzkm); 
                
                obj.covariances(:,:,idx) = obj.covariances(:,:,idx) ...
                                           - Pxzk*Kgain' - Kgain*Pxzk' ...
                                           + Kgain*Pzzk*Kgain'; 
            end
        end

        function [cond_mean, cond_cov] = getCondProps(obj)
            % For a given state of the filter, this function calculates and
            % outputs the conditional mean and covariance of the gaussian
            % mixture

             cond_mean = zeros(obj.stateSize,1); 
             
             for idx = 1:obj.numWeights
                cond_mean = cond_mean + obj.means(:,:,idx) * obj.weights(idx); 
             end

             cond_cov = zeros(size(obj.covariances(:,:,1)));
             
             for idx = 1:obj.numWeights
                cond_cov = cond_cov + obj.weights(idx)*...
                (obj.covariances(:,:,idx) + (obj.means(:,:,idx) - cond_mean)*(obj.means(:,:,idx) - cond_mean)');
             end
        end

        function out = distribution(obj)
            % For each term of the state, this function generates a 2000
            % term PDF distribution across an interval that covers the
            % -4sigma and +4sigma from the smallest to the largest means,
            % enabling the plotting of GM distributions. 

            minMeanCov.mean = zeros(1, obj.stateSize); 
            minMeanCov.cov  = zeros(1, obj.stateSize);
            
            maxMeanCov.mean = zeros(1, obj.stateSize); 
            maxMeanCov.cov  = zeros(1, obj.stateSize); 


            for idx = 1:obj.stateSize
                minMeanCov.mean(idx)    = min(obj.means(idx, :));
                associatedCov           = diag(obj.covariances(:,:,find(obj.means(idx,:) == minMeanCov.mean(idx), 1))); 
                minMeanCov.cov(idx)     = associatedCov(idx);

                maxMeanCov.mean(idx)    = max(obj.means(idx, :));
                associatedCov           = diag(obj.covariances(:,:,find(obj.means(idx,:) == maxMeanCov.mean(idx), 1))); 
                maxMeanCov.cov(idx)     = associatedCov(idx);
            end
            
            numElements= 2000; 
            xEval = zeros(obj.stateSize, numElements);

            for idx = 1:obj.stateSize
                xEval(idx, :) = linspace(minMeanCov.mean(idx) - 4*sqrt(minMeanCov.cov(idx)), ...
                    maxMeanCov.mean(idx) + 4*sqrt(maxMeanCov.cov(idx)),numElements); 
            end

            dist = zeros(obj.stateSize,numElements); 
            
 
            for jdx = 1:obj.stateSize
               for idx = 1:obj.numWeights          
                    covars = diag(obj.covariances(:,:,idx));
                    dist(jdx,:) = dist(jdx,:) + obj.weights(idx) * normpdf(xEval(jdx,:), obj.means(jdx, 1, idx), sqrt(covars(jdx)));
               end
            end

            out.distributions = dist; 
            out.xEval = xEval; 
        end

        function randomSample = sampleDistribution(obj)
            mu = rand(1);
            csum = cumsum(obj.weights);
            
            jdx = find(mu<=csum, 1,"first"); 

            randomSample = obj.means(:,1,jdx) + randn(1, obj.stateSize) *  diag(chol(obj.covariances(:,:,jdx))');
        end
    end
end