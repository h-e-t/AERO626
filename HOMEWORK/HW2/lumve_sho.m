    function [Xhat_hist, P_hist, K_hist, innov_hist, S_hist] = lumve_sho(T, Z, R, x0, P0, omega, Q)
    % Sequential LUMVE (BLUE/Kalman) for SHO with position-only measurements.
    % Inputs:
    %   T      : 1xN time vector (monotonic)
    %   Z      : 1xN position measurements
    %   R      : scalar measurement variance (or 1x1)
    %   x0     : 2x1 initial state [x; xdot]
    %   P0     : 2x2 initial covariance
    %   omega  : natural frequency (default 1)
    %   Q      : 2x2 process-noise cov per step (default zeros)
    % Outputs:
    %   Xhat_hist : 2xN posterior state estimates
    %   P_hist    : 2x2xN posterior covariances
    %   K_hist    : 2x1xN gains
    %   innov_hist: 1xN innovations z - H x_b
    %   S_hist    : 1xN innovation variances
    
    if nargin < 6 || isempty(omega), omega = 1; end
    if nargin < 7 || isempty(Q), Q = zeros(2); end
    
    N  = numel(T);
    H  = [1 0];
    I2 = eye(2);
    
    % Preallocate
    Xhat_hist  = zeros(2, N);
    P_hist     = zeros(2, 2, N);
    K_hist     = zeros(2, 1, N);
    innov_hist = zeros(1, N);
    S_hist     = zeros(1, N);
    
    % Initialization
    x_p   = x0(:);
    P_p   = P0;
    
    % Helper: STM(dt)
    Phi_of = @(dt) [ cos(omega*dt),           (1/omega)*sin(omega*dt);
                    -omega*sin(omega*dt),     cos(omega*dt)           ];
    
    for k = 1:N
        % --- Time step ---
        if k == 1
            dt = T(1);          % assume t0 = 0; adjust if needed
        else
            dt = T(k) - T(k-1);
        end
        Phi = Phi_of(dt);
    
        % --- Predict ---
        x_b = Phi * x_p;
        P_b = Phi * P_p * Phi.' + Q;   % include Q if you model process noise
    
        % --- Update (LUMVE/BLUE) ---
        % Innovation and covariance
        y    = Z(k) - H * x_b;                 % innovation
        S    = H * P_b * H' + R;              % scalar here
        K    = (P_b * H') / S;                % 2x1 gain
    
        % Posterior (use Joseph form for numerical robustness)
        x_hat = x_b + K * y;
        % P_hat = (I2 - K*H) * P_b * (I2 - K*H).' + K*R*K.';
        P_hat = (I2 - K*H)*P_b;
    
        % --- Save ---
        Xhat_hist(:, k) = x_hat;
        P_hist(:, :, k) = P_hat;
        K_hist(:, :, k) = K;
        innov_hist(k)   = y;
        S_hist(k)       = S;
    
        % --- Roll forward ---
        x_p = x_hat;
        P_p = P_hat;
    end
    end