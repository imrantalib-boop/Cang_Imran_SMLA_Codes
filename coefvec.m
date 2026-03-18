function C=coefvec(f,N,L)
syms t
W=weight(L);
Psi=funvec(N,L);
KW=((f)*(Psi)*(W)) % Product
for i=0:N
     H=ort(i)
    aa=(KW(i+1))
    aa_1=matlabFunction(aa);
 C(i+1)=((1)/(H))*integral(aa_1,0,1)
end
    C=double(C)
end



















% Y = KW(:,2)
% q = KW(:,3)
% K=  KW(:,4)
% L= KW(:,5)
%c= KW(:,6)

%     for i=0:M
%      H=ort(M)
%      aa=(KW(i+1))
%  C(i+1)=1/H*int(aa,0,1)
%     end
% %end
% r = KW(:,1);
% x =  [0.1 0.2 0.3];
%B=trapz(x,r)


