function [E1,E2]=Esym(n,W,G1,G2,b,const,V2,m1,m2)
im=sqrt(-1);
e1= (n+0.5)*const*b/m1-im*G1;
e2= -(n-0.5)*const*b/m2+W-im*G2;
a=sqrt((e1-e2)^2+8*n*V2*b);
E1=0.5*(e1+e2+a);
E2=0.5*(e1+e2-a);
end