clear 
clc
% 
% load("firstSave.mat")
% 
% syms I t tau dt sw sf sbw sbf real
% 
% assume(dt, 'real')
% 


% Z33 = zeros(3,3); 
% I3  = eye(3); 
% I5  = eye(5); 
% 
% F = [0 0 0 -1 0;
%      0 0 0 0 -1;
%      0 1 0 0 0;
%      0 0 0 0 0;
%      0 0 0 0 0];



% STM = I5 + F *tau + 1/2 * F^2. * tau^2;
% 
% 
% Qc   = diag([sw sf 0 sbw sbf]);
% 
% expr = STM * Qc * STM';
% 
% Qd   = int(expr, tau, 0, dt);

% function Qd = calc_Qd(obj, sw, sf, sbw, sbf, dt)
%     Z33 = zeros(3,3); 
% 
%     sw = sw ^ 2; 
%     sf = sf ^ 2; 
%     sbw = sbw ^ 2; 
%     sbf = sbf ^ 2; 
% 
%     Qd = ... 
%     [(sbw*dt^3)/3 + sw*dt,                         Z33,                        Z33, -(dt^2*sbw)/2,           Z33;
%                       Z33,       (sbf*dt^3)/3 + sf*dt,  (sbf*dt^4)/8 + (sf*dt^2)/2,           Z33, -(dt^2*sbf)/2;
%                       Z33, (sbf*dt^4)/8 + (sf*dt^2)/2, (sbf*dt^5)/20 + (sf*dt^3)/3,           Z33, -(dt^3*sbf)/6;
%            -(dt^2*sbw)/2,                         Z33,                         Z33,        dt*sbw,           Z33;
%                       Z33,              -(dt^2*sbf)/2,               -(dt^3*sbf)/6,           Z33,        dt*sbf];
% end

function magMeas = measurement(attitude, defaultMag)
    arguments
        attitude 
        defaultMag = [27.5550, -2.4169, -16.0849]'
    end

    magMeas = quat2dcm(attitude)*defaultMag; 
    
end

defaultMag = [27.5550, -2.4169, -16.0849]'

yaw = 90
rotation = eul2quat(deg2rad([yaw 0 0]), "ZYX");

meas = measurement(rotation)



quiver(0,0,defaultMag(1),defaultMag(2), DisplayName="0 yaw")
hold on
quiver(0,0,meas(1),meas(2), DisplayName=yaw+" yaw")
xlim([-30, 30])
ylim([-30, 30])
legend
axis("equal")