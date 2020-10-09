function G = hpfiltereasy(Y,l)
T=length(Y);
% For simplicity, assume T is at least 5
F=zeros(T);

F(1,1)=l+1;
F(1,2)=-2*l;
F(1,3)=l;
F(2,1)=-2*l;
F(2,2)=5*l+1;
F(2,3)=-4*l;
F(2,4)=l;

for i=3:T-2
    F(i,i-2)=l;
    F(i,i-1)=-4*l;
    F(i,i)=6*l+1;
    F(i,i+1)=-4*l;
    F(i,i+2)=l;
end

F(T-1,T-3)=l;
F(T-1,T-2)=-4*l;
F(T-1,T-1)=5*l+1;
F(T-1,T)=-2*l;
F(T,T-2)=l;
F(T,T-1)=-2*l;
F(T,T)=l+1;

G=F\Y;

