% 30 equally spaced data points and various N are considered 
tic
clc
clear all
close all

syms t
N=8;%Neurons
L = 1;
v = t^9; % Analytical solution
t_values = linspace(0, 1, 30); % Training data
% Ordfer of the differential equations
    gamma2 = 2; 
    gamma = 5/2;
    gamma1 = 1/2;
    
    u = 4*t^9 - (196608/(36465*sqrt(pi)))*t^(17/2) + 72*t^7 + (49152/(143*sqrt(pi)))*t^(13/2); % Source term
    phi = funvec(N, L); % Chebyshev Polynomial generation
     U = double(right_hand_matrix(t_values, u)); % Source term at training points
     E = double(matrix_descretize_points(L, N, gamma2, gamma, gamma1, t_values)); % Discretized matrix

     % Compute optimal weights
     MPGI1 = double(inv(transpose(E)*E)*transpose(E)); % Moore-Penrose inverse
     weight = transpose(MPGI1 * U);

%%%%%%%%%%%%%%%%%%%%%%Absolute Error%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
NNS= weight*phi; %Neural Network solution

tn=0:0.1:1
 Error=subs(abs(v-NNS),symvar(t),tn);
 fprintf('Values in scientific notation:\n');
for i = 1:length(Error)
    fprintf('% .2e\n', Error(i));  % .4e means scientific notation with 4 decimal digits
end
 Maxerror=vpa(max(Error))

toc



