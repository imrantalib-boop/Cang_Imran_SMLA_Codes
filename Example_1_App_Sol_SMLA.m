% Solution computed using SMLA algorithm considering N=25 and 30 equally
% spaced data points
clc
clear all
close all

syms t
L = 1;
v = t^9; % Analytical solution
t_values = linspace(0, 1, 30); % Training data
% Ordfer of the differential equations
    gamma2 = 2; 
    gamma = 5/2;
    gamma1 = 1/2;
    
    u = 4*t^9 - (196608/(36465*sqrt(pi)))*t^(17/2) + 72*t^7 + (49152/(143*sqrt(pi)))*t^(13/2); % Source term
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Computational of Total Loss%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=25;
phi = funvec(N, L); % Chebyshev Polynomial generation
U = double(right_hand_matrix(t_values, u)); % Source term at training points
E = double(matrix_descretize_points(L, N, gamma2, gamma, gamma1, t_values)); % Discretized matrix
MPGI1 = double(inv(transpose(E)*E)*transpose(E)); % Moore-Penrose inverse
%%%%%%%%%%%%%%%%%%%%%%Neural Network Solution%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 weight = transpose(MPGI1 * U)
NNS= weight*phi; %Neural Network solution

hold on

fplot(v,[0,1],'LineWidth',2)

tn_vals = linspace(0,1,20);             
NNS_fun = matlabFunction(NNS);         
plot(tn_vals, NNS_fun(tn_vals),'o', ...
     'LineStyle','none','MarkerSize',7,'LineWidth',1.5)

hold off

xlabel('$t$','interpreter','latex','fontsize',20)
ylabel('$v(t)$','Interpreter','latex','FontSize',20)

legend({'Analytical solution','SMLA, $N=25$'},...
       'Location','best','FontSize',20,'Interpreter','latex')

box on
set(gca,'FontSize',20)




