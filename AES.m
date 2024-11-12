load('C1209.mat');
h=figure;
subplot(211);plot(Time,V, 'LineWidth',2); grid on
xlabel('Time(hr)');
ylabel('Meas. Voltage(V)');
subplot(212);plot(Time,I, 'LineWidth',2); grid on
xlabel('Time(hr)');
ylabel('Meas. Current(A)');
sgtitle('Batt:C1209');

% **********************Question 1**********************
dt = [diff(Time); 0]; % Time intervals
charge_current = I(I > 0); % Charging current
discharge_current = I(I < 0); % Discharging current
Q_c = sum(charge_current .* dt(I > 0)); % Charge capacity in Ah
Q_d = abs(sum(discharge_current .* dt(I < 0))); % Discharge capacity in Ah

% Display results
disp('******Question 1******');
fprintf('Charge Capacity (Qc): %.4f Ah\n', Q_c);
fprintf('Discharge Capacity (Qd): %.4f Ah\n', Q_d);


% **********************Question 2**********************
SOC = ones(length(Time), 1); % Initialize SOC array with initial SOC=1

% Compute SOC over time
for k = 2:length(Time)
    if I(k) > 0
        Q = Q_c; % Charging capacity
    else
        Q = Q_d; % Discharging capacity
    end
    % Calculate SOC at each time step using dt
    SOC(k) = SOC(k-1) + (dt(k) * I(k)) / (3600 * Q);
end

% Plot Question 2 SOC vs. Time
figure;
plot(Time, SOC, 'LineWidth', 2);
xlabel('Time (hr)');
ylabel('State of Charge (SOC)');
title('Question 2 SOC vs. Time');
grid on;


% **********************Question 3**********************
% Define epsilon for scaling
epsilon = 0.175;
s_scaled = epsilon + (1 - 2 * epsilon) * SOC; % Apply linear scaling

% Define the matrix P for both models
% Combined model
P_combined = [ones(size(s_scaled)), 1 ./ s_scaled, s_scaled, log(s_scaled), log(1 - s_scaled)];

% Combined+3 model
P_combined3 = [ones(size(s_scaled)), 1 ./ s_scaled, 1 ./ (s_scaled.^2), 1 ./ (s_scaled.^3), ...
                    1 ./ (s_scaled.^4), s_scaled, log(s_scaled), log(1 - s_scaled)];

% Calculate the least squares solution for both models
% Combined model
c_combined = (P_combined' * P_combined) \ (P_combined' * V);
% c_combined = round(c_combined, 4); % Round to 4 decimal places

% Combined+3 model
c_combined3 = (P_combined3' * P_combined3) \ (P_combined3' * V);
% c_combined3 = round(c_combined3, 4); % Round to 4 decimal places

% Display the parameters
disp('Combined model parameters:');
disp(c_combined);

disp('Combined+3 model parameters:');
disp(c_combined3);


% **********************Question 4**********************
% Number of data points
N = length(V);

% Model parameters count
M_combined = length(c_combined);
M_combined3 = length(c_combined3);

% Predicted voltage for Combined model
V_combined_hat = P_combined * c_combined;

% Calculate Error Metrics for Combined model
V_mean = mean(V);
R2_combined = (1 - norm(V - V_combined_hat)^2 / norm(V - V_mean)^2) * 100;
MaxError_combined = max(abs(V - V_combined_hat));
RMS_combined = norm(V - V_combined_hat) / sqrt(N - M_combined);

% Display Results for Combined model
fprintf('--- Combined Model ---\n');
fprintf('R² Fit: %.2f%%\n', R2_combined);
fprintf('Max Error: %.4f V\n', MaxError_combined);
fprintf('Root-Mean Square Error (RMS): %.4f V\n', RMS_combined);

% Predicted voltage for Combined+3 model
V_combined3_hat = P_combined3 * c_combined3;

% Calculate Error Metrics for Combined+3 model
R2_combined3 = (1 - norm(V - V_combined3_hat)^2 / norm(V - V_mean)^2) * 100;
MaxError_combined3 = max(abs(V - V_combined3_hat));
RMS_combined3 = norm(V - V_combined3_hat) / sqrt(N - M_combined3);

% Display Results for Combined+3 model
fprintf('\n--- Combined+3 Model ---\n');
fprintf('R² Fit: %.2f%%\n', R2_combined3);
fprintf('Max Error: %.4f V\n', MaxError_combined3);
fprintf('Root-Mean Square Error (RMS): %.4f V\n', RMS_combined3);

% **********************Question 5**********************
% Define SOC range from 0 to 1 with a step of 0.01
s = (0:0.01:1);
s = epsilon + (1 - 2 * epsilon) * s; 
% Define the parameter array (replace with your actual estimated values)
k = c_combined3; % replace with your array values

% Combined+3 Model Calculation using the array values
V_combined3 = k(1) + k(2)./s + k(3)./s.^2 + k(4)./s.^3 + k(5)./s.^4 + ...
       k(6) * s + k(7) * log(s) + k(8) * log(1 - s);

% Plot the Combined+3 Model
figure;
plot(s, V_combined3, 'r-', 'LineWidth', 2);

% Add labels and title
xlabel('State of Charge (SOC)');
ylabel('Open Circuit Voltage (V_{o})');
title('OCV-SOC Curve for Combined+3 Model');
grid on;
