%approximate solution is computed for example_1 using spectral tau approach for N=25 
clear all
close all
clc
syms t f c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 d0 d1 d2 d3 d4 d5 d6 d6 d8 d9
% Coefficient of differential equations
a1=1; a2=1; a3=1;
L=1;
% Scale level of the numerical scheme
N=25;
% order of the derivatives
gamma2=2; gamma=5/2; gamma1=1/2;
V=t^9;
%Source terms
u=4*t^9-(196608/(36465*sqrt(pi)))*t^(17/2)+72*t^7+(49152/(143*sqrt(pi)))*t^(13/2);
% Initial conditions
v0=subs(V,t,0)
v1=subs(diff(V,1),t,0)
v2=subs(diff(V,2),t,0)


PHI= funvec(N,L);% Legendre function vector
G= coefvec(u,N,L);% Coefficients of the known source term

% Fractional derivative approximation using  Operational matrices of derivatives in Caputo sense

Dgamma2= ope_derr(L,gamma2,N); %Order 2
Dalpha1= ope_derr(L,1,N); % Order 1
Dgamma= ope_derr(L,gamma,N); % Order 3/2
Dgamma1= ope_derr(L,gamma1,N);% Order 1/2

%Residual calculation
C=[c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 d0 d1 d2 d3 d4 d5]%N=25
 R2= (C*Dgamma2+C*Dgamma-2*C*Dgamma1+4*C-G)*PHI; 
 
 % Computation of m-n equations using Tau method
 
 Innerproduct=int(R2*PHI,0,1);
 
 % Computation of equations using initial conditions
 A1=C*subs(PHI,t,0)-v0;
 A2=C*Dalpha1*subs(PHI,t,0)-v1;
 A3=C*Dgamma2*subs(PHI,t,0)-v2

 % calculation of coefficients for solving system of algebraic equations
 [c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 d0 d1 d2 d3 d4 d5]=solve(Innerproduct(1),Innerproduct(2),Innerproduct(3),Innerproduct(4),Innerproduct(5),Innerproduct(6),Innerproduct(7),Innerproduct(8),Innerproduct(9),Innerproduct(10),Innerproduct(11),Innerproduct(12),Innerproduct(13),Innerproduct(14),Innerproduct(15),Innerproduct(16),Innerproduct(17),Innerproduct(18),Innerproduct(19),Innerproduct(20),Innerproduct(21),Innerproduct(22),Innerproduct(23),A1,A2,A3)
 Coef=[c0 c1 c2 c3 c4 c5 c6 c7 c8 c9 b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 d0 d1 d2 d3 d4 d5]
 Vapp=Coef*PHI
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
hold on
fplot(V,[0,1])
fplot(Vapp,[0,1])
  hold off
  xlabel('$t$','interpreter','latex','fontsize',25)
  ylabel('$v(t)$','Interpreter','latex','FontSize',25)
 legend({'Analytical solution','Spectral Tau Method, $N=25$'},'location','best','fontsize',25,'Interpreter','latex')
 box on
set(gca,'FontSize',20)
 
 
 
 
 