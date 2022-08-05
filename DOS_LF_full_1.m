function [f1,f2]=DOS_LF_full_1(omega,B,m1,m2,const,G1,G2,W,Gap,Ev)
%Following function calculates the DOS using both the solutions separately
%There was a sign issue in DOS_LF_full which has been rectified here
im=sqrt(-1);
mu0=W*m2/(m1+m2);
omega=omega-mu0;
const=1.157*10^-1;%this is e\hbar/m_e  in units of meV*m_e/T
if G1==0 
    eta=im*W/20000;
    f1=1;f2=1;
else
    f1=0;f2=0;
binv=B^-1;
Eb=const*B/m1;
Eb2=const*B/m2;
mm=m1*m2/(m2-m1);
mp=m1*m2/(m2+m1);
Ebm=const*B/mm;
Ebp=const*B/mp;
Ehyb=Eb*Ev;
%contribution from the symmetric combination
a=-Eb*Eb2;
b=Eb*(Eb2+W-im*G2)+im*Eb2*G1-(mu0+omega)*Ebm-2*Ehyb;
c=-im*G1*(W+Eb2-im*G2)-(mu0+omega)*(Eb2+W-im*G1-im*G2)+(mu0+omega)^2+Ehyb;
d1=-Ebm;d1=d1/a;
d2=2*(omega+mu0)-Eb2-W+im*(G1+G2);d2=d2/a;

x1=-b/(2*a);
x2=sqrt(b^2-4*a*c)/(2*a);x2=x2*sign(imag(x2));
xp=x1+x2;
xm=x1-x2;
 
c1=d1/2;
c2=(d1*(xp+xm)/2+d2)/(xm-xp);
cp=c1-c2;
cm=c1+c2;

% a=mu0/Eb+((m1-m2)*omega+im*(m1*G1-m2*G2))/(2*const*B);
% b=sqrt(((m1+m2)*omega+im*(m1*G1+m2*G2))^2-m1*m2*Gap^2)/(2*const*B);
% b=b*sign(imag(b));
% xp=a+b;
% xm=a-b;
% b=sqrt(((m1+m2)*omega+im*(m1*G1+m2*G2))^2-m1*m2*Gap^2);
% b=b*sign(imag(b));
Pfact=-B;
f1=imag(Pfact*im*pi*(cp*sign(imag(xp))*exp(2*pi*im*sign(imag(xp))*xp)));
f2=imag(Pfact*im*pi*(sign(imag(xm))*cm*exp(2*pi*im*sign(imag(xm))*xm)));

%contribution from the asymmetric combination
a=-Eb*Eb2;
b=Eb*(W-im*G2)+Eb2*(im*G1+Eb)-(mu0+omega)*Ebm-2*Ehyb;
c=-(Eb+im*G1)*(W-im*G2)-(mu0+omega)*(W-im*G1-im*G2-Eb)+(mu0+omega)^2+Ehyb;
d1=-Ebm;d1=d1/a;
d2=2*(omega+mu0)+Eb-W+im*(G1+G2);d2=d2/a;

x1=-b/(2*a);
x2=sqrt(b^2-4*a*c)/(2*a);x2=x2*sign(imag(x2));
xp=x1+x2;
xm=x1-x2;
 
c1=d1/2;
c2=(d1*(xp+xm)/2+d2)/(xm-xp);
cp=c1-c2;
cm=c1+c2;

% a=mu0/Eb+((m1-m2)*omega+im*(m1*G1-m2*G2))/(2*const*B);
% b=sqrt(((m1+m2)*omega+im*(m1*G1+m2*G2))^2-m1*m2*Gap^2)/(2*const*B);
% b=b*sign(imag(b));
% xp=a+b;
% xm=a-b;
% b=sqrt(((m1+m2)*omega+im*(m1*G1+m2*G2))^2-m1*m2*Gap^2);
% b=b*sign(imag(b));
Pfact=-B;
f1=f1+imag(Pfact*im*pi*(cp*sign(imag(xp))*exp(2*pi*im*sign(imag(xp))*xp)));
f2=f2+imag(Pfact*im*pi*(sign(imag(xm))*cm*exp(2*pi*im*sign(imag(xm))*xm)));
end
end