classdef BParticleFilter < handle
    % This class implements the Bootstrap Particle Filter 

    properties
        means       (:,1,:)
        covariances (:,:,:)
        weights     (1,:)
        stateSize   (1,1)       
        numWeights  (1,1)


    end

    methods
        function obj = BParticleFilter(means,covariances,numParticles)            
            arguments
                means 
                covariances 
                numParticles   
            end
            
            obj.means           = means; 
            obj.covariances     = covariances; 
            obj.weights         = numParticles; 
            obj.stateSize       = length(means(:,1)); 

        end


    end

end