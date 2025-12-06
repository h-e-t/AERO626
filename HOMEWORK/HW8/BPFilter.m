classdef BPFilter < handle

    properties
       numParticles
       
       ParticleSwarm
       weights  
       
       resamplingThreshold
       resamplingEnabled = false

       stateSize

       resampleCount = 0
   end

    methods
        function obj = BPFilter(numParticles, sampleDist, resamplingEnabled, resamplingThreshold)
            arguments
                numParticles 
                sampleDist 
                resamplingEnabled    = false
                resamplingThreshold  = 0 
            end
            
            
            obj.numParticles = numParticles; 
            
            obj.weights = ones(1,numParticles)/numParticles;
            
            obj.ParticleSwarm = zeros(1,numParticles);
            for idx = 1:numParticles
                obj.ParticleSwarm(idx) = sampleDist.sampleDistribution(); 
            end

            if resamplingEnabled
                assert(resamplingThreshold > 0)

                obj.resamplingEnabled = true;
                obj.resamplingThreshold = resamplingThreshold; 
            end
        end

        function propagate(obj, F, Pww)
            % F is passed as non linear function of the state with noise
            % included so emulate the process correctly 

            for idx = 1:obj.numParticles
                obj.ParticleSwarm(idx) = F(obj.ParticleSwarm(idx), Pww); 
            end
        end

        function update(obj, H, Pvv, zk)
            for idx = 1:obj.numParticles
                Zxk = H(obj.ParticleSwarm(idx)); 
                obj.weights(idx) = normpdf(zk, Zxk, Pvv) * obj.weights(idx); 
            end

            obj.weights = obj.weights / sum(obj.weights);
            
            if obj.resamplingEnabled && obj.calculateEffectiveParticles() < obj.resamplingThreshold 
                obj.resample();
            end
        end

        function resample(obj)
            csum = cumsum(obj.weights); 
            
            for idx = 1:obj.numParticles
                zeta = rand(); 

                jdx = find(zeta <= csum, 1, "first");
                
                obj.ParticleSwarm(idx) = obj.ParticleSwarm(jdx); 
            end

            obj.weights = ones(1, obj.numParticles) * 1/obj.numParticles; 

            obj.resampleCount = obj.resampleCount + 1; 
        end
        
        function [mean, cov] = getProps(obj)
            mean = sum(obj.ParticleSwarm .* obj.weights); 

            cov = 0;
            for idx = 1:obj.numParticles
                cov = cov + obj.weights(idx)*(obj.ParticleSwarm(idx) - mean)^2; 
            end
        end

        function Neff = calculateEffectiveParticles(obj)
            Neff = 1 / sum(obj.weights.^2);
        end
   end
end