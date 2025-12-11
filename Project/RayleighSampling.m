function [x_meas,v_meas] = RayleighSampling(position, velocity)
    
    % Sigma defined by 2.5 50% CEP by NEO6M gps
    r_sigma = 1.177; 
   
    r = r_sigma * sqrt(-2*log(1-rand())); 

    theta = 2*pi*rand(); 

    dr = [r * cos(theta); r * sin(theta); 1.5*r_sigma*randn()] ;
    dv = [.1*randn(); .1*randn(); .2*randn()]; 
    
    x_meas = position + dr; 
    v_meas = velocity + dv; 
end