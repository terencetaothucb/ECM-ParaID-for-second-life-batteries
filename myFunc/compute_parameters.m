function [R1, R2, C, estimated_Q_mean, mape] = compute_parameters(theta, ocv_vector, soc_vector, Q_ref)
    % Compute physical parameters and estimate battery capacity with MAPE evaluation
    % Inputs:
    %   theta - Parameter vector from parameter identification
    %   ocv_vector - Open Circuit Voltage (OCV) vector
    %   soc_vector - State of Charge (SOC) vector
    %   Q_ref - The reference capacity value for mape calculation
    % Outputs:
    %   R1 - Estimated resistance R1
    %   R2 - Estimated resistance R2
    %   C - Estimated capacitance C
    %   estimated_Q_mean - Mean estimated capacity Q
    %   mape - Mean Absolute Percentage Error (MAPE) for capacity estimation

    % Compute resistance R1 directly from theta
    R1 = theta(1);

    % Compute capacitance C using derived formula
    C = 1 / (theta(3) / theta(4) + theta(2) + theta(1) * theta(4));

    % Compute resistance R2 using C and theta
    R2 = 1 / (theta(4) * C);

    % Compute differences
    docv = diff(ocv_vector); % Compute differences in OCV
    dsoc = diff(soc_vector); % Compute differences in SOC
    beta_1 = docv ./ dsoc; % Ratio of differences as an alternative gradient
    beta_1(isnan(beta_1)) = 1e-8; % Handle potential NaN values due to division by zero
    estimated_Q = - beta_1 * theta(4) / theta(3); % Recompute estimated Q with beta_1

    % Compute the mean of estimated capacities, ignoring invalid values
    estimated_Q_mean = mean(estimated_Q(isfinite(estimated_Q))); % Exclude NaN and Inf values

    % Handle edge case where estimated_Q_mean is NaN
    if isnan(estimated_Q_mean)
        estimated_Q_mean = 1e6; % Assign a large default value if estimation fails
    end

    % Compute Mean Absolute Percentage Error (MAPE) for capacity estimation
    mape = abs((estimated_Q_mean - Q_ref) / Q_ref) * 100; % Percentage error relative to true Q
end