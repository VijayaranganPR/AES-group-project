clear;
close all;
clc;

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
% Identify the indices where the current I is positive (charging) and negative (discharging)
charging_indices = find(I > 0);
discharging_indices = find(I < 0);

% Calculate the total charging and discharging time
charging_time = Time(charging_indices(end)) - Time(charging_indices(1));
discharging_time = Time(discharging_indices(end)) - Time(discharging_indices(1));

% Compute the average absolute current during charging and discharging
average_charging_current = mean(abs(I(charging_indices(1):charging_indices(end))));
average_discharging_current = mean(abs(I(discharging_indices(1):discharging_indices(end))));

% Calculate the total charge for charging and discharging phases
Qc = charging_time * average_charging_current;
Qd = discharging_time * average_discharging_current;

% Display results
disp('**********************Question 1**********************');
fprintf('Charge Capacity (Qc): %.4f Ah\n', Qc);
fprintf('Discharge Capacity (Qd): %.4f Ah\n', Qd);

% **********************Question 2**********************
SOC = zeros(length(I), 1);
SOC(1) = 1;

for i = 2:length(I)
    delta = Time(i) - Time(i-1);
    if(I(i))<0
        SOC(i) = SOC(i-1) + (delta * I(i))/Qd;
    else
        SOC(i) = SOC(i-1) + (delta * I(i))/Qc;
    end

    if (SOC(i) < 0)
        SOC(i) = 0;
    elseif(SOC(i) > 1)
        SOC(i) = 1;
    end
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
P_combined = [ones(size(s_scaled)), 1 ./ s_scaled, s_scaled, log(s_scaled), log(1 - s_scaled), I];

% Combined+3 model
P_combined3 = [ones(size(s_scaled)), 1 ./ s_scaled, 1 ./ (s_scaled.^2), 1 ./ (s_scaled.^3), ...
                    1 ./ (s_scaled.^4), s_scaled, log(s_scaled), log(1 - s_scaled), I];

% Calculate the least squares solution for both models
% Combined model
c_combined = (P_combined' * P_combined) \ (P_combined' * V);
c_combined = round(c_combined, 4); % Round to 4 decimal places

% Combined+3 model
c_combined3 = (P_combined3' * P_combined3) \ (P_combined3' * V);
c_combined3 = round(c_combined3, 4); % Round to 4 decimal places

disp('**********************Question 3**********************')
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
disp("**********************Question 4**********************")
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
s_raw = (0:0.01:1);
s = epsilon + (1 - 2 * epsilon) * s_raw; 
k = c_combined3;

% Combined+3 Model Calculation using the array values
V_combined3 = k(1) + k(2)./s + k(3)./s.^2 + k(4)./s.^3 + k(5)./s.^4 + ...
       k(6) * s + k(7) * log(s) + k(8) * log(1 - s);

% Plot the Combined+3 Model
figure; hold on; box on;
plot(s_raw, V_combined3, 'r-', 'LineWidth', 2);


% Add labels and title
xlabel('State of Charge (SOC)');
ylabel('Open Circuit Voltage (V_{o})');
title('Question 5 OCV-SOC Curve for Combined+3 Model');
grid on;

% **********************Question 6**********************
fprintf('\n*********Question 6a*********\n')
SOC_levels = [0.25, 0.50, 0.75];
OCV_values = interp1(s_raw, V_combined3, SOC_levels, 'linear'); 
fprintf('OCV when the battery SOC is 25%%: %.4f V\n', OCV_values(1));
fprintf('OCV when the battery SOC is 50%%: %.4f V\n', OCV_values(2));
fprintf('OCV when the battery SOC is 75%%: %.4f V\n', OCV_values(3));

fprintf('*********Question 6b*********\n')
% Target OCV value
target_OCV = 4.0;
% Interpolating SOC at 4V
SOC_at_4V = interp1(V_combined3, s_raw, target_OCV, 'linear');

fprintf('SOC at 4V (zero current): %.4f\n', SOC_at_4V);

fprintf('*********Question 6c*********\n')
% Target OCV value with a discharge current of 1A
Vo_OCV_3_8V = 3.8;

Reff = c_combined3(end);
discharge_current6c = -1;
target_OCV_3_8V = Vo_OCV_3_8V - discharge_current6c * Reff;
% Interpolating SOC at 3.8V
SOC_at_3_8V = interp1(V_combined3, s_raw, target_OCV_3_8V, 'linear');

fprintf('SOC at %.4fV with 1A discharge: %.4f\n', target_OCV_3_8V, SOC_at_3_8V);



% **********************Question 7**********************
% from question 4
fprintf('\n*********Question 7a*********\n')

% Check for duplicates in the SOC array
[SOC_unique, idx_unique] = unique(SOC, 'stable'); % Get unique SOC values and their indices
V_combined3_hat_unique = V_combined3_hat(idx_unique); % Corresponding OCV values without duplicates

OCV_values_7a = interp1(SOC_unique, V_combined3_hat_unique, SOC_levels, 'linear');
fprintf('OCV when the battery SOC is 25%%: %.4f V\n', OCV_values_7a(1));
fprintf('OCV when the battery SOC is 50%%: %.4f V\n', OCV_values_7a(2));
fprintf('OCV when the battery SOC is 75%%: %.4f V\n', OCV_values_7a(3));

fprintf('*********Question 7b*********\n')
% Target OCV value
target_OCV = 4.0;
% Interpolating SOC at 4V
SOC_at_4V = interp1(V_combined3_hat_unique, SOC_unique, target_OCV, 'linear');

fprintf('SOC at 4V (zero current): %.4f\n', SOC_at_4V);

fprintf('*********Question 7c*********\n')
% Target OCV value with a discharge current of 1A
Vo_OCV_3_8V = 3.8;

Reff = c_combined3(end);
discharge_current6c = -1;
target_OCV_3_8V = Vo_OCV_3_8V - discharge_current6c * Reff;
% Interpolating SOC at 3.8V
SOC_at_3_8V = interp1(V_combined3_hat_unique, SOC_unique, target_OCV_3_8V, 'linear');

fprintf('SOC at %.4fV with 1A discharge: %.4f\n', target_OCV_3_8V, SOC_at_3_8V);

