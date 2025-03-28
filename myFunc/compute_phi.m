function [phi1, phi2, phi3, phi4, phi5] = compute_phi(lambda_0, lambda_1, time, current, V_deviation, dt)
    % Compute state-space filter outputs for the parameter identification model
    % Inputs:
    %   lambda_0, lambda_1 - Filter parameters
    %   time - Time vector
    %   current - Current signal vector I(t)
    %   V_deviation - Voltage deviation signal Ṽ(t)
    %   dt - Time step vector
    % Outputs:
    %   phi1, phi2, phi3, phi4, phi5 - Filtered output signals

    % Define the state-space matrices
    A = [0, 1; -lambda_0, -lambda_1]; % State matrix
    B = [0; 1];                       % Input matrix

    % Define output and direct transmission matrices for each filter
    C_phi1 = [-lambda_0^2, -lambda_1 * lambda_0]; 
    D_phi1 = lambda_0;
    C_phi2 = [0, lambda_0];                   
    D_phi2 = 0;
    C_phi3 = [lambda_0, 0];                      
    D_phi3 = 0;
    C_phi4 = [0, -lambda_0];                    
    D_phi4 = 0;
    C_phi5 = [-lambda_0, 0];                   
    D_phi5 = 0;

    % Initialize filtered output vectors
    phi1 = zeros(length(time), 1); 
    phi2 = zeros(length(time), 1); 
    phi3 = zeros(length(time), 1); 
    phi4 = zeros(length(time), 1); 
    phi5 = zeros(length(time), 1); 

    % Initialize state vectors
    x1 = [0; 0]; 
    x2 = [0; 0]; 
    x3 = [0; 0]; 
    x4 = [0; 0]; 
    x5 = [0; 0]; 

    % Iterate over each time step to compute filter outputs
    for k = 1:length(time)-1
        % Update states for filters related to current signal
        x1_dot = A * x1 + B * current(k); 
        x2_dot = A * x2 + B * current(k);
        x3_dot = A * x3 + B * current(k);

        % Compute filtered outputs
        phi1(k) = C_phi1 * x1 + D_phi1 * current(k);
        phi2(k) = C_phi2 * x2 + D_phi2 * current(k);
        phi3(k) = C_phi3 * x3 + D_phi3 * current(k);

        % Update states using Euler integration
        x1 = x1 + x1_dot * dt(k);
        x2 = x2 + x2_dot * dt(k);
        x3 = x3 + x3_dot * dt(k);

        % Update states for filters related to voltage deviation signal
        x4_dot = A * x4 + B * V_deviation(k);
        x5_dot = A * x5 + B * V_deviation(k);

        % Compute filtered outputs
        phi4(k) = C_phi4 * x4 + D_phi4 * V_deviation(k);
        phi5(k) = C_phi5 * x5 + D_phi5 * V_deviation(k);

        % Update states using Euler integration
        x4 = x4 + x4_dot * dt(k);
        x5 = x5 + x5_dot * dt(k);
    end

    % Compute final values of the filtered outputs
    phi1(end) = C_phi1 * x1 + D_phi1 * current(end);
    phi2(end) = C_phi2 * x2 + D_phi2 * current(end);
    phi3(end) = C_phi3 * x3 + D_phi3 * current(end);
    phi4(end) = C_phi4 * x4 + D_phi4 * V_deviation(end);
    phi5(end) = C_phi5 * x5 + D_phi5 * V_deviation(end);
end