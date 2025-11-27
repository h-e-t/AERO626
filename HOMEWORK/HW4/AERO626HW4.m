%% Problem 2
clc 
clear all

rand(626); 

mx_0    = [1;0];           % m, m/s
Pxx_0   = diag([4, 1]);   %(m/s)^2
Hx       = [1 0];

Pvv = .1^2; % m^2

function m = STM(t)
    m = [cos(t) sin(t);
          -sin(t) cos(t)];
end
true_time = 0:.1:20; 
true_meas = [];
for idx = true_time
    true_meas = [true_meas STM(idx) * [1;0]];  
end




t = 1:1:20;
zk = [];
for idx = t
    zk = [zk STM(idx) * [1;0]];  
end


priori_mean             = zeros(2,length(t)); 
posterior_mean          = zeros(2,length(t));

priori_covariance       = zeros(2,length(t)); 
posterior_covariance    = zeros(2,length(t));
innovation_covariance   = zeros(1,length(t));

error_history       = [];
arbitrary_t         = [0];
state_cov_history   = [sqrt(diag(Pxx_0))];
meas_cov_history    = [];

state_history       = [];

mx_0    = [1;0];          % m, m/s
Pxx_0   = diag([4, 1]);   %(m/s)^2
Hx      = [1 0];
Hv      = 1;

Pvv = .1^2; % m^2

mx_m    = mx_0; 
pxx_m   = Pxx_0; 


for idx = 1:length(t)
    Fx = STM(1); 
    
    % Propogate 
    mx_m        = Fx*mx_m;
    pxx_m       = Fx*pxx_m*Fx'; % no process noise

    % Update 
    mz_m        = Hx * mx_m; 

    pxz_m       = pxx_m * Hx'; 
    pzz_m       = Hx*pxx_m*Hx' + Hv*Pvv*Hv'; 

    K           = pxz_m / pzz_m; 

    mx_p        = mx_m + K*(zk(idx) - mz_m);
    pxx_p       = pxx_m - pxz_m*K' - K*pxz_m' + K*pzz_m*K';
     
    % RECORD 

    priori_mean(:, idx) = mx_m;
    posterior_mean(:, idx) = mx_p;
    
    priori_covariance(:, idx)       = sqrt(diag(pxx_m));
    posterior_covariance(:, idx)    = sqrt(diag(pxx_p));

    innovation_covariance(idx)    = sqrt(pzz_m); 
    
    arbitrary_t         = [arbitrary_t t(idx) t(idx)];
    state_cov_history   = [state_cov_history priori_covariance(:,idx) 
                            posterior_covariance(:, idx)];
    

    % SETUP NEXT ITERATION
    mx_m    = mx_p;
    pxx_m   = pxx_p;
end

%%

figure(1)
tiledlayout(2,1)
nexttile
plot(arbitrary_t,  state_cov_history(1,:), "Color","red")
title("Position Covariance")

hold on
plot(arbitrary_t, -state_cov_history(1,:),"Color","red")
ylabel("\bf{\pm 1 \sigma (m)} ", Interpreter="tex")

nexttile
plot(arbitrary_t, state_cov_history(2,:), "Color","blue")
title("Velocity Covariance")
hold on
plot(arbitrary_t, -state_cov_history(2,:), "Color","blue")
xlabel("Time (S)")
ylabel("\bf{\pm 1 \sigma (m/s)}", Interpreter="tex")

figure

plot(t,ones(1,length(t))*sqrt(Pvv), Color='red', DisplayName="Pvv \pm1\sigma" )
hold on
plot(t,-ones(1,length(t))*sqrt(Pvv), Color='red', HandleVisibility='off')


plot(t, innovation_covariance, Color="magenta", DisplayName="Pzz \pm1\sigma" )
plot(t, -innovation_covariance, Color="magenta", HandleVisibility='off')

ylabel("m^2")
xlabel("Time (S)")


legend()
