function estimate = LUMVEEstimator(z, t, F, Ht, STM, W, x0, w0)
    arguments
        z (:,:) double % Measurements List  1xm
        t (1,:) double % time reference     1xm
        F (:,:) double % Model Matrix       
        Ht (:,:) double % Measurement Mapping Matrix
        STM function_handle % State Transition Matrix Function 
        W  (:,:) double = NaN % Weighting matrix 
        x0 (:,:) double = NaN % Previous state
        w0 (:,:) double = NaN % Previous state covariance
    end

    H = [];
    
    for idx = 1:length(t)
        stm = STM(t(idx)); 
        
        Hk = Ht * stm; 
        
        H = [H; Hk]; 
    end

    

    if ~any(isnan(W)) && ~any(isnan(x0)) 
        weights = diag(W);
        estimate = (H'*weights*H + w0)\(H'*weights*z + w0*x0); 
    else if ~any(isnan(W))
        weights = diag(W);
        estimate = (H'*weights*H)\(H'*weights*z);
    elseif ~any(isnan(x0))  
        estimate = (H'*H + w0)\(H'*z + w0*x0); 
    else
        estimate = (H'*H)\(H'*z);
    end
end