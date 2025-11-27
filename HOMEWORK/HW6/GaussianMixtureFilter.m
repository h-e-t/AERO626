classdef GaussianMixtureFilter < handle
    % This class implements the linear Gaussian mixture filter for an 
    % arbitrary state size and arbitrary state dynamics

    properties
        means       (:,1,:)
        covariances (:,:,:)
        weights     (1,:)
        stateSize   (1,1)       
        numWeights  (1,1)        
    end

    methods
        function obj = GaussianMixtureFilter(means,covariances,weights)
            obj.means           = means; 
            obj.covariances     = covariances; 
            obj.weights         = weights; 
            obj.stateSize       = length(means(:,1)); 
            obj.numWeights      = length(weights); 
        end

        function obj = predictPDF(obj, F, Pww)
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

        function obj = updatePDF(obj, H, Pvv, Z)
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
                    maxMeanCov.mean(idx) + 4*sqrt(maxMeanCov.cov(idx)) ,numElements); 
            end

            dist = zeros(obj.stateSize,numElements); 

            for jdx = 1:obj.stateSize
               for idx = 1:length(obj.weights)
                    covars = diag(obj.covariances(:,:,idx));
                    
                    dist(jdx,:) = dist(jdx,:) + obj.weights(idx) * normpdf(xEval(jdx,:), obj.means(jdx, idx), sqrt(covars(jdx)));
               end
            end

            out.distributions = dist; 
            out.xEval = xEval; 
        end


        function randomState = sampleDistribution(obj)
            dist = obj.distribution(); 

            samples = dist.distribution; 

            randi(length(sample(1,:)));

            randomState = zeros(obj.stateSize, 1);
            
            for i = 1:obj.stateSize
                
            end
        end
    end
end