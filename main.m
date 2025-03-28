clear
%% Data Loading
% User input for the dataset name
dataset_name = 'Cylind21'; 
Q_nominal = 2.1; % Nominal capacity in Ah
% Specify SOC and SOH
selected_SOC = 25; % Specified SOC value
%%
selected_SOH_col = 1; % Specified SOH column index
% Selected SOH column index: 1, 8, 65 
% to reproduce the reported results
%%

% Construct filenames based on the dataset name
raw_data = ['data_', dataset_name, '_timeserials.xlsx'];
OCV_SOC_data = ['OCV_SOC_', dataset_name, '_Interp.xlsx'];

% Load the sheet names of the raw data file
sheet_names = sheetnames(raw_data);

% Load the first sheet to calculate the number of columns
data_sample = readtable(raw_data, 'Sheet', sheet_names{1}, 'VariableNamingRule', 'preserve');
num_columns = width(data_sample);

% Compute the maximum SOH index
col_width = num_columns - 2; 

% Convert numeric SOC to string format
selected_SOC_name = ['SOC', num2str(selected_SOC)]; 

% Verify and load the selected SOC sheet
if ~ismember(selected_SOC_name, sheet_names)
    error(['Invalid SOC: ', selected_SOC_name, '. Ensure it matches a sheet name in the file.']);
end
data = readtable(raw_data, 'Sheet', selected_SOC_name, 'VariableNamingRule', 'preserve');

% Extract SOH value from the column header
SOH_name = data.Properties.VariableNames{selected_SOH_col + 2};
tokens = regexp(SOH_name, 'SOH=([\d.]+)', 'tokens');
if isempty(tokens) || isempty(tokens{1})
    error('Unable to extract SOH value from the column header.');
end
SOH = str2double(tokens{1}{1});

% Extract all SOH values
SOH_values = zeros(1, col_width);
for k = 1:col_width
    SOH_column_name = data_sample.Properties.VariableNames{k + 2};
    tokens = regexp(SOH_column_name, 'SOH=([\d.]+)', 'tokens');
    if ~isempty(tokens) && ~isempty(tokens{1})
        SOH_values(k) = str2double(tokens{1}{1});
    else
        error(['Unable to extract SOH value from column header: ', SOH_column_name]);
    end
end

% Find the maximum SOH value and its index
[max_SOH_value, Index_SOH_max] = max(SOH_values);

% Compute the true capacity
Q_true = SOH * Q_nominal;

% Extract time, current, and voltage data
time = data{:, 1};        % Time
current = data{:, 2};     % Current
voltage = data{:, selected_SOH_col + 2}; % Voltage
dt = diff(time);          % Time intervals

%% Load Interpolated OCV Data
interpolated_OCV_SOC_data = readtable(OCV_SOC_data);
interpolated_OCV_SOC_data = interpolated_OCV_SOC_data(interpolated_OCV_SOC_data.BatteryIndex == Index_SOH_max, :);

% Create interpolation functions
soc2ocv = @(soc_input) interp1(interpolated_OCV_SOC_data.SOC_Interp, interpolated_OCV_SOC_data.OCV_Interp, soc_input);
ocv2soc = @(ocv_input) interp1(interpolated_OCV_SOC_data.OCV_Interp, interpolated_OCV_SOC_data.SOC_Interp, ocv_input);

% Initialize SOC and OCV vectors
init_soc = ocv2soc(voltage(1));
soc_vector = init_soc * ones(length(voltage), 1);
ocv_vector = voltage(1) * ones(length(voltage), 1);

% Compute SOC and OCV for each timestep
soc_k = init_soc;
for k = 1:length(voltage)-1
    soc_k = soc_k + dt(k) / 3600 / Q_nominal * current(k) * 100; % SOC update based on current
    soc_vector(k+1) = soc_k;
    ocv_vector(k+1) = soc2ocv(soc_k);
end

%% Optimization Loop
lambda0_list = linspace(1e-6, 1e1, 300); % Range for λ0
lambda1_list = linspace(1e-6, 1e1, 300); % Range for λ1

% Initialize structure for storing the best parameters
best_params = struct('lambda_0', 0, 'lambda_1', 0, 'mape', inf);
feer = zeros(length(lambda0_list),length(lambda1_list));
% Grid search over λ0 and λ1
for i = 1:length(lambda0_list)
    for j = 1:length(lambda1_list)
        lambda_0 = lambda0_list(i);
        lambda_1 = lambda1_list(j);

        % Compute filters and parameters
        [z, phi_matrix] = compute_filters(lambda_0, lambda_1, time, voltage, current, ocv_vector, dt);
        theta = phi_matrix \ z; % Least-squares estimate of θ

        % Compute physical parameters
        [R1, R2, C, estimated_Q, mape] = compute_parameters(theta, ocv_vector, soc_vector, Q_nominal);
        
        feer(i, j) = mape;

        % Update the best parameters if MAPE improves
        if mape < best_params.mape
            best_params.lambda_0 = lambda_0;
            best_params.lambda_1 = lambda_1;
            best_params.mape = mape;
            best_params.R1 = R1;
            best_params.R2 = R2;
            best_params.C = C;
            best_params.estimated_Q = estimated_Q;
        end
    end
end

%% Display Selected Inputs
disp('--- Selected Inputs ---');
disp(['Selected SOC value: ', num2str(selected_SOC), '%']);
disp(['Selected SOH column index: ', num2str(selected_SOH_col)]);
disp(['Selected SOH value to test: ', num2str(SOH)]);
disp(['OCV-SOC curve comes from BOL SOH: ', num2str(max_SOH_value)]);

%% Display Optimization Results
disp('--- Optimization Results ---');
disp(['Best λ0: ', num2str(best_params.lambda_0)]);
disp(['Best λ1: ', num2str(best_params.lambda_1)]);
disp(['R1: ', num2str(best_params.R1), ' Ohms']);
disp(['R2: ', num2str(best_params.R2), ' Ohms']);
disp(['C: ', num2str(best_params.C), ' Farads']);
disp(['Tau: ', num2str(best_params.R2 * best_params.C), ' s']);
disp(['Q_true: ', num2str(Q_true), ' Ah']);
disp(['Q_nominal: ', num2str(Q_nominal), ' Ah']);
disp(['Estimated Q: ', num2str(best_params.estimated_Q), ' Ah']);
disp(['MAPE to Q_nominal: ', num2str(best_params.mape), ' %']);
disp(['MAPE to Q_true: ', num2str(abs(best_params.estimated_Q - Q_true) / Q_true * 100), ' %']);