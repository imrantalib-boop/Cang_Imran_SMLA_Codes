function F = right_hand_matrix(t_values,f)
syms t
F=double(subs(f,t_values));

% Reshape F to be a column matrix
F = reshape(F, [], 1);

% Add the initial conditions z(0) = 0 and z(1) = exp(-1)
    F(end + 1, 1) = 0;              % z(0) = 0
    F(end + 1, 1) = 0;             % z'(0) = 0
    F(end + 1, 1) = 0;             % z''(0)=0
    F;
end