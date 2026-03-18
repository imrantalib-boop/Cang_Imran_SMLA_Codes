%Example1 approximate solution using operational matrices approach at N=25
%at various N
clear all
close all
clc
tic
syms t
L=1;
gamma=5/2;
gamma1=1/2
gamma2=2;
v=t^9;
G=4*t^9-(196608/(36465*sqrt(pi)))*t^(17/2)+72*t^7+(49152/(143*sqrt(pi)))*t^(13/2);
v0=subs(v,t,0);
v1=subs(diff(v,t),t,0);
v2=subs(diff(v,t,2),t,0)
f1=v0+v1*t;
N=25;
    Psi=funvec(N,L)
   
H=double(ope_Int(gamma2,N,L))
Q1=double(ope_derr(L,gamma,N));
Q2=double(ope_derr(L,gamma1,N))
B=double(H*Q1-2*H*Q2+4*H);
F1=zeros(1,N+1)
F2=coefvec(G,N,L)
D= double(F1*Q1-2*F1*Q2+4*F1-F2);   
 A=1;
 LL=lyap(A,B,D)
 Uapp=LL*H*Psi+ F1*Psi;%solution computed using OMA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 hold on
fplot(v,[0,1])
fplot(Uapp,[0,1])
  hold off
  xlabel('$t$','interpreter','latex','fontsize',32)
  ylabel('$v(t)$','Interpreter','latex','FontSize',32)
 legend({'Analytical solution','Operational matrices approach, $N=25$'},'location','best','fontsize',35,'Interpreter','latex')
 box on
set(gca,'FontSize',32)
 