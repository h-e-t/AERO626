%% Problem 1
clc
clear
load("data_HW03.mat")

k1 = 2.5;
k2 = 3.7;
m  = 1.5;
h  = 5.4;


w_n = sqrt((k1 + k2)/m) ;

F = [0 1; -w_n^2 0];

function stm = STM(dt, w)
    arguments
        dt 
        w = 2.033060090930254
    end
    stm = [cos(w*(dt))    1/w*sin(w*(dt)); 
           -w*sin(w*(dt)) cos(w*(dt))]; 

end

function Ht = Jac(xs, h)
    arguments
        xs 
        h = 5.4
    end
    x = xs(1);
    v = xs(2);
    Ht = [x/sqrt(x^2 + h^2), 0;
          v/sqrt(x^2 + h^2) - x^2*v/(sqrt(x^2 + h^2))^3, x/sqrt(x^2+h^2)];
end

function meas = measurementMapping(xs, h)
    arguments
        xs 
        h = 5.4
    end
    x = xs(1);
    v = xs(2);
    meas = [sqrt(x^2 + h^2);
        x*v/sqrt(x^2 + h^2)]; 
end

xs0 = [4;.2];

x0 = [3;0];

dbx0 = [0 ; 0];
Pbxx0 = diag([1000 100]);
Pvvi = eye(2); 
Zp1 = [Range; RangeRate]; 


f1 = figure; 
tiledlayout(4,2)

t0 = Time(1); 

colors = ["rx", "bx", "gx", "mx"];
for n = 1:4

    lam = Pbxx0\dbx0; 
    LAM = inv(Pbxx0); 

    residuals = [];

    for l = 1:length(Time)
        phi = STM(Time(l) - t0);
        xsl = phi * xs0; 
        Ht = Jac(xsl); 
        h =  measurementMapping(xsl);
        zl = Zp1(:,l); 

        residuals = [residuals zl - h];

        lam = lam + (Ht*phi)'/Pvvi * (zl - h);
        LAM = LAM + (Ht*phi)'/Pvvi * (Ht * phi); 
        

    end

    dhx0 = LAM\lam;

    xs0  = xs0 + dhx0;

    dbx0 = dbx0 - dhx0;
    nexttile

    plot(residuals(1,:), colors(n), 'DisplayName',"Iteration "+n, LineWidth=2); 
    if n == 1 
        title("Range Residual")
    end

    nexttile

    plot(residuals(2,:), colors(n), 'DisplayName',"Iteration "+n, LineWidth=2); 
    legend()

    if n == 1 
        title("Range Rate Residual")
    end
    


 
end
fprintf("State after 4 iterations = (%.6d, %.6d)^T", xs0)



%% Problem 2
clear
clc
close all
load("data_HW03.mat")

function d = fallingBodyDeriv(~, s)
    rho0    = 3.4e-3; % lb/ft3
    kp      = 22000;  % ft
    g       = 32.2;   % ft/s2.

    x = s(1);
    v = s(2);
    beta = s(3);

    d = [v;1/2 * rho0 * exp(-x/kp) * v^2/beta-g; 0];
end

function F = dynJac(s)
    rho0    = 3.4e-3; % lb/ft3
    kp      = 22000;  % ft
    
    x = s(1);
    v = s(2);
    beta = s(3);

    t21 = -rho0/(2*kp) * exp(-x/kp)*v^2/beta;
    t22 = rho0*exp(-x/kp)*v/beta;
    t23 = -rho0/2*exp(-x/kp)*v^2/beta^2; 

    F = [0    1    0;
         t21, t22, t23;
         0    0    0]; 
end

function Ht = mJac(x)
    r1      = 1000;   % ft
    r2      = 10;     % ft

    Ht = [(x(1) + r2)/sqrt(r1^2 + (x(1) + r2)^2); 0; 0];

end

function [t, xs] = fallingBodyTraj(times, x0)
    [t, xs] = ode45(@fallingBodyDeriv, times, x0(:)); 
end

function [t,xs] = fbTraj_Stm(times,x0)
    function d = fbstm(t,s)
        x_dot   = fallingBodyDeriv(t,s); 
        stm_dot = dynJac(s) * reshape(s(end-8:end), [3,3]);

        d = [x_dot; stm_dot(:)];
    end

    [t, xs] = ode45(@fbstm, times, x0(:));
end

function m = measMap(x)
    r1      = 1000;   % ft
    r2      = 10;     % ft

    m = sqrt(r1^2 + (x(1) + r2)^2);
end



Pvv     = 100;    % ft2



xs0 = [ 105000; 
        -5700; 
        1750]; 
% 
xs0 = [ 80000; 
        -6500; 
        2250]; 

I3 = eye(3); 

[t, xs] = fbTraj_Stm(T, [xs0;I3(:)]);

dbx0 = 0;

figure; 
numResiduals = 8;
tiledlayout(numResiduals,2)

maxIter = 1000; 
tol = .2; 
RMS = [];
correctionMagnitude = []; 
estimatedStates = [xs0]; 
tileIDX = 1;  

for n = 1:maxIter

    [t, xs] = fbTraj_Stm(T, [xs0;I3(:)]);

    PHI = reshape(xs(:,4:12)', 3, 3, []);
    XS = xs(:,1:3); 

    lam = 0; 
    LAM = 0; 

    residuals = [];

    for l = 1:length(T)
        phi = PHI(:,:,l);
        xsl = XS(l,:); 
        Ht  = mJac(xsl)'; 
        h   = measMap(xsl);
        zl  = Z(l); 

        residuals = [residuals zl - h];

        lam = lam + (Ht*phi)'/Pvv * (zl - h);
        LAM = LAM + (Ht*phi)'/Pvv * (Ht * phi); 


    end

    dhx0 = LAM\lam;

    xs0  = xs0 + dhx0;

    dbx0 = dbx0 - dhx0;

    RMS = [RMS sqrt(1/length(T) * sum(residuals.^2/Pvv))]; 
    correctionMagnitude = [correctionMagnitude norm(dhx0)];
    estimatedStates = [estimatedStates xs0];

    disp("Iter " + n + " processed")

    nexttile(tileIDX, [1,1])
    tileIDX = tileIDX + 2; 
    plot(residuals(1,:), 'x', 'DisplayName',"Iteration "+n, LineWidth=2); 
    legend

    if norm(dhx0) < tol
        break;
    end
end

nexttile([numResiduals,1])
plot(T, Z,'bx', DisplayName="Measurements")
hold
idx = 1; 
for s = estimatedStates
    [t, TRAJ] = fallingBodyTraj(T, s);
    plot(t, TRAJ(:,1), DisplayName="Trajectory Estimate "+idx)
    idx = idx + 1; 
end
legend

saveas(gcf, "residualntraj3.pdf")


f2= figure;
tiledlayout(1,2)
nexttile
plot(1:length(RMS), RMS)
ylabel(gca, "RMS", Interpreter="latex")
xlabel(gca, "Iteration", Interpreter="latex")

nexttile 
plot(1:length(RMS),correctionMagnitude)
xticks(1:length(RMS))
ylabel(gca, "||$\delta \hat{x_0}$||", Interpreter="latex")
xlabel(gca, "Iteration", Interpreter="latex")

saveas(f2, "RMS and correction3.pdf")


