
format long
clc
clear all
close all

syms t
L = 1;
v = t^9; % Analytical solution
t_values = linspace(0, 1, 30); % Training data
% Ordfer of the derivatives
    gamma2 = 2; 
    gamma = 5/2;
    gamma1 = 1/2;
    
    u = 4*t^9 - (196608/(36465*sqrt(pi)))*t^(17/2) + 72*t^7 + (49152/(143*sqrt(pi)))*t^(13/2); % Source term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Computational of Total Loss%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize arrays to store results
 m_values = 5:5:25; % Neuron counts to test
 losses = zeros(size(m_values)); % Store total losses

 for i = 1:length(m_values)
     m = m_values(i);
 
     phi = funvec(m, L); % Chebyshev Polynomial generation
     U = double(right_hand_matrix(t_values, u)); % Source term at training points
     E = double(matrix_descretize_points(L, m, gamma2, gamma, gamma1, t_values)); % Discretized matrix

    % Compute optimal weights
    MPGI1 = double(inv(transpose(E)*E)*transpose(E)); % Moore-Penrose inverse
     weight = transpose(MPGI1 * U); % Optimal weights
     L_function = abs(E * transpose(weight) - U); 
    total_loss = sum(L_function); % Total loss
 
     losses(i) = total_loss;

    fprintf('m = %d, Total Loss = %.4e\n', m, total_loss);
 end

% Plotting
% figure;
plot(m_values, losses, '-o', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('$N$', 'Interpreter', 'latex', 'FontSize', 32);
ylabel('Total Loss', 'Interpreter', 'latex', 'FontSize', 32);
legend('Loss vs. Neuron count', 'Interpreter', 'latex', 'FontSize', 36);
box on;
set(gca, 'FontSize', 32);
set(gca,'YScale','log');