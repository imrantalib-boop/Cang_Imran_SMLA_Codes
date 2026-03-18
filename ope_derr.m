%chebyshev derivative operational matrix
function D=ope_derr(L,v,N)
ceil_neu=ceil(v);
D=zeros(N+1);
for i=ceil_neu:N
    Svji=zeros(1,N+1);
    for j=0:N
                          if   j==0
                              etaj=2;
                          else
                              etaj=1;
                          end
    for k=ceil_neu:i
        
n1=(-1)^(i-k)*2*i*factorial(i+k-1)*gamma(k-v+0.5);
d1=etaj*L^(v)*gamma(k+0.5)*factorial(i-k)*gamma(k-v-j+1)*gamma(k+j-v+1);
Svji(k-ceil_neu+1)=n1/d1;    
    end
    D(i+1,j+1)=sum(Svji);
    end
end

D
end