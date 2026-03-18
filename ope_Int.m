function H=ope_Int(sigma,N,L)
syms k l ccc
ccc= zeros(1,N+1);
H=zeros(N+1);
for i=0:N
    for j=0:N
                          if   j==0
                              etaj=2;
                          else
                              etaj=1;
                          end

                          for k=1:i
    aa=((-1)^(i-k))*2*i*L^(sigma)*factorial(i+k-1)*gamma(k+sigma+1/2);
    bb= etaj*gamma(k+1/2)*factorial(i-k)*gamma(k+sigma-j+1)*gamma(k+j+sigma+1);
   % aa=((-1)^(i-k))*2*factorial(i)*factorial(k)*L^(sigma)*factorial(i+k-1)*gamma(k+sigma+1/2);
   % bb= etaj*gamma(k+1/2)*factorial(i-k)*factorial(2*k)*gamma(k+sigma-j+1)*gamma(k+j+sigma+1)
   ccc(k+1)=(aa/bb)
                          end
                          H(i+1,j+1)=sum(ccc)     

    end


end

H
end


