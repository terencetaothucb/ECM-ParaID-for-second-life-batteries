function [z, phi_matrix] = compute_filters(lambda_0, lambda_1, time, voltage, current, ocv_vector, dt)
    % Compute the filtered signal z(t) and the phi matrix for parameter estimation
    % Inputs:
    %   lambda_0, lambda_1 - Filter parameters
    %   time - Time vector
    %   voltage - Measured voltage vector V(t)
    %   current - Measured current vector I(t)
    %   ocv_vector - Open Circuit Voltage (OCV) vector
    %   dt - Time step vector
    % Outputs:
    %   z - Filtered signal z(t)
    %   phi_matrix - Matrix of filtered signals \(\phi_1\) to \(\phi_5\)

    % Define state-space matrices for the z(t) filter
    A_z = [0, 1; -lambda_0, -lambda_1]; % State matrix
    B_z = [0; 1];                       % Input matrix
    C_z = [-lambda_0^2, -lambda_1 * lambda_0]; % Output matrix
    D_z = lambda_0;                     % Direct transmission matrix

    % Initialize state variables
    x_z = [0; 0]; % Initial state vector for z(t)
    z = zeros(length(voltage), 1); % Initialize z output vector
    V_deviation = zeros(length(voltage), 1); % Initialize voltage deviation vector

    ocv_current = ocv_vector(1); % Initial OCV value

    % Compute filtered signal z(t)
    for k = 1:length(voltage)-1
        % Calculate voltage deviation as difference between measured voltage and OCV
        V_deviation(k) = voltage(k) - ocv_current;

        % State-space model: compute state derivative
        x_dot_z = A_z * x_z + B_z * V_deviation(k);

        % Compute filtered signal z(t)
        z(k) = C_z * x_z + D_z * V_deviation(k);

        % Update state using Euler integration
        x_z = x_z + x_dot_z * dt(k);

        % Update current OCV value
        ocv_current = ocv_vector(k);
    end

    % Compute final values for the last time step
    V_deviation(end) = voltage(end) - ocv_current;
    z(end) = C_z * x_z + D_z * V_deviation(end);

    % Compute \(\phi_1\) to \(\phi_5\) using a helper function
    [phi1, phi2, phi3, phi4, phi5] = compute_phi(lambda_0, lambda_1, time, current, V_deviation, dt);

    % Construct the \(\phi\) matrix by combining filtered signals
    phi_matrix = [phi1, phi2, phi3, phi4, phi5];
end