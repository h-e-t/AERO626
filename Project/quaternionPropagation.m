attitude = quaternion(1,0,0,0); 

sigma = deg2rad([0,0,5])';
dt    = 1; 



function new_q = quatPropagation(q, w, dt)

    sigma = w * dt;
    q     = q.compact'; 

    sigma = deg2rad([0,0,5]');

    sigma_norm = norm(sigma);

    sin_sig_over_2 = sin(sigma_norm/2); 
    cos_sig_over_2 = cos(sigma_norm/2); 
    

    exp_OMEGA_sigma = [cos_sig_over_2,         sin_sig_over_2/2*sigma';
                       sin_sig_over_2/2*sigma, cos_sig_over_2*eye(3) - sin_sig_over_2/2*skewSym(sigma)];
    
    new_q = (exp_OMEGA_sigma * q)'; 
    

end

quatPropagation(attitude, sigma, 5)
