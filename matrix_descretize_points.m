function A= matrix_descretize_points(L,m,alpha,alpha1,alpha2,t_values)%m is nunmber of neurons and alpha(s) is order of operational matrix
syms t
DD2 =ope_der(L,alpha,m);% alpha is highest order of derivative and m is number of neurons in hidden layers L is interval end point[0,L].
DD52 =ope_der(L,alpha1,m);
DD12 =ope_der(L,alpha2,m);
phi=funvec(m,L);

% Loop over Chebyshev polynomials (i) and interval points (j)
 for i = 0:m
     for j = 1:length(t_values)
         t = t_values(j);

% Compute the matrix entry A(i,j)

A(j, i+1)= subs(DD2(i+1),t)+subs(DD52(i+1),t)-2*subs(DD12(i+1),t)+4*subs(phi(i+1),t);

     end
 end
 % Calculate Legendre polynomials at t = 0 using the shift_Chebyshev function
     Chebyshev_at_t0 = funvec(m,L)
     values_at_t0 = subs(Chebyshev_at_t0, 0);
     values_at_t1 = subs(ope_der(L,1,m), 0);
     values_at_t2 = subs(ope_der(L,2,m), 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
     A(end+1, :) = values_at_t0;
     A(end+1, :) = values_at_t1;
     A(end+1, :) = values_at_t2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 A=double(A);
end