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

%%
clc
% Define the initial vector
% firstVector = [1; 0; 0];

% Create a DCM for a 45-degree rotation of a vector about the Z-axis
g = [0 0 9.81]';
dcm45degZ = angle2dcm(0, deg2rad(5), deg2rad(0), "XYZ")'; 
% skewSym(dcm45degZ * g) * [] 


% Convert DCM to quaternion
% expectedQuaternion = dcm2quat(dcm45degZ)

% Convert quaternion back to DCM
% expectedQuaternionDCM = quat2dcm(expectedQuaternion);

% Transform the initial vector
% deg45Vector = expectedQuaternionDCM * firstVector

% Display the transformed vector
% disp(deg45Vector);






%%
clc
clear

% tic 
attitude = quaternion(1,0,.5,0); 


% quat2rotm(attitude)
% quat2dcm(attitude)
rotQuat = quaternion(quatexp(.5*[0,0,0,deg2rad(45)]));
% % 
% for i = 1:400
% rotQuat = quaternion(quatexp(.5*[0,0,deg2rad(5),deg2rad(5)]));
% rotQuat = quaternion(quatexp(.5*[0,0,0,deg2rad(45)]));
smallRot = rotQuat  ;
% end 
% % disp("Operation took: "+toc+" seconds")
quat2rotm(smallRot)*[1;0;0]





%%
a = [1 2 3]';
b = [3 2 1]';

skewSym(a)*b
-skewSym(b)*a
