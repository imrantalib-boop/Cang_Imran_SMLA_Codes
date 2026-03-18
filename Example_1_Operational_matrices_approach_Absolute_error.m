%Example1 using operational matrices approach; absolute error is computed
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
N=6;
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
 tn=0:0.1:1
 Error=vpa(subs(abs(v-Uapp),symvar(t),tn));
 fprintf('Values in scientific notation:\n');
for i = 1:length(Error)
    fprintf('% .2e\n', Error(i));  
end
 Maxerror=max(Error)

 toc
 
 