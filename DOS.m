function f=DOS(omega,n_c,b,m1,m2,const,G1,G2,V2,W)
im=sqrt(-1);
if G1==0 
    eta=im*W/10000;
else
    eta=0;
end

f=0;
    binv=b^-1;
    es0=0.5*const*b/m1-im*G1;ea0=-0.5*const*b/m2+W-im*G2;
    f=f-1/pi*imag(1./(omega+eta-es0))-1/pi*imag(1./(omega+eta-ea0));
    %Following sum ensures correct result for omega>W
    for k=1:n_c
        [E1,E2]=Esym(k,W,G1,G2,b,const,V2,m1,m2);
        f=f-1/pi*imag(1./(omega+eta-E1))-1/pi*imag(1./(omega+eta-E2));
        [E1,E2]=Easym(k,W,G1,G2,b,const,V2,m1,m2);
        f=f-1/pi*imag(1./(omega+eta-E1))-1/pi*imag(1./(omega+eta-E2));
    end
   
    %following sum is for trickier omegas lying below W
    if omega<W &&  floor((W-omega)/(const/(m2*binv)))>n_c
       nf=floor((W-omega)/(const/(m2*binv)));
       a=f;c=1;
       [E1,E2]=Esym(nf,W,G1,G2,b,const,V2,m1,m2);
       f=f-1/pi*imag(1./(omega+eta-E1))-1/pi*imag(1./(omega+eta-E2));
       [E1,E2]=Easym(nf,W,G1,G2,b,const,V2,m1,m2);
       f=f-1/pi*imag(1./(omega+eta-E1))-1/pi*imag(1./(omega+eta-E2));
       while abs(a-f)>10^-4
           a=f;
           [E1,E2]=Esym(nf+c,W,G1,G2,b,const,V2,m1,m2);
           f=f-1/pi*imag(1./(omega+eta-E1))-1/pi*imag(1./(omega+eta-E2));
           [E1,E2]=Easym(nf+c,W,G1,G2,b,const,V2,m1,m2);
           f=f-1/pi*imag(1./(omega+eta-E1))-1/pi*imag(1./(omega+eta-E2));
       if nf-c>n_c
           [E1,E2]=Esym(nf-c,W,G1,G2,b,const,V2,m1,m2);
           f=f-1/pi*imag(1./(omega+eta-E1))-1/pi*imag(1./(omega+eta-E2));
           [E1,E2]=Easym(nf-c,W,G1,G2,b,const,V2,m1,m2);
           f=f-1/pi*imag(1./(omega+eta-E1))-1/pi*imag(1./(omega+eta-E2));
       end
       c=c+1;
       end
   end
    f=f*b;
end
