%testing the k-sum in Kubo formula


clc;clear all;
im=sqrt(-1);t1=100;
%following way of defining the paramters are from FW paper but scaled down
Gap=2;%in meV
W=27*Gap;%in meV
m1=1;m2=10*m1;mplus=m2*m1/(m2+m1);mminus=m2*m1/(m2-m1);%in m_e
mu0=W*m2/(m1+m2);%in meV
kconst=6.06*10^-3;
kf=kconst*sqrt(2*mplus*W);%in units of a0^-1
p1=4;P=linspace(0.3,.9,p1);
v=Gap/(2*kf);%in units of meV.a0
Omega=-100:0.01:100;


for p2=1:p1
G1=P(p2)*Gap;G2=0.1*Gap;
Gm=G1+G2;gm=G1-G2;
Econst=13605;%hbar^2/(2.me.a0^2)
ek1=@(k) Econst*k.^2/(m1)-mu0;ek2=@(k) Econst*k.^2/(m2)-W+mu0;
delta=@(k) v*(k);
Ek1=@(k) 0.5*(ek1(k)-ek2(k)-im*Gm);
Ek2=@(k) 0.5*sqrt((ek1(k)+ek2(k)-im*gm).^2+4*abs(delta(k).^2));

Ekp=@(k) real(Ek1(k)+Ek2(k));taup=@(k) abs(-1./(2*imag(Ek1(k)+Ek2(k))));
Ekm=@(k) real(Ek1(k)-Ek2(k));taum=@(k) abs(-1./(2*imag(Ek1(k)-Ek2(k))));

vx1=@(k) 0.5*Econst*(2*k/m1-2*k/m2);
vx2=@(k) real(0.25*(2*(ek1(k)+ek2(k)-im*gm).*(Econst*2*k/m1+Econst*2*k/m2)+8*v^2*k)./sqrt((ek1(k)+ek2(k)-im*gm).^2+8*abs(delta(k)).^2));
vxp=@(k) vx1(k)+vx2(k);vxm=@(k) vx1(k)-vx2(k);
dosp=@(k,omega) -1/pi*imag(1./(omega-Ek1(k)-Ek2(k)));
dosm=@(k,omega)-1/pi*imag(1./(omega-Ek1(k)+Ek2(k)));

tht_cont=2/3;

%FOLLOWING FUNCTIONS ARE USED TO DETERMINE THE BOUNDARY OF PARTITIONED 2D INTEGRAL 
Ep=@(k) ek1(k)-ek2(k)+abs(ek1(k)+ek2(k));
Em=@(k) ek1(k)-ek2(k)-abs(ek1(k)+ek2(k));
Fp=@(k,omega) tht_cont*k.*dosp(k,omega).^2.*vxp(k).^2;
Fm=@(k,omega) tht_cont*k.*dosm(k,omega).^2.*vxm(k).^2;
%THIS LOOP DEFINES THE k-integration as a function of omega
for j=1:max(size(Omega))
fp=@(k) Fp(k,Omega(j));
Ip(j)=integral(fp,0,1);
fm=@(k) Fm(k,Omega(j));
Im(j)=integral(fm,0,1);
end
%t=linspace(0.01,.01,t1);t(1)=0;
t=logspace(-3,1,t1);t(1)=0;
C=zeros(1,t1);
for i=1:t1
    if t(i)==0
    C(i)=Kubo_0T(G1);
else
    beta=t(i)^-1;
    Om_c=9/beta;
    dfp=@(omega) beta*exp(-beta*(omega))./(1+exp(-beta*(omega))).^2;
    dfm=@(omega) beta*exp(beta*(omega))./(1+exp(beta*(omega))).^2;
    Fp1=@(omega) interp1(Omega,Ip,omega).*dfp(omega);
    Fm1=@(omega) interp1(Omega,Im,omega).*dfm(omega);
    C(i)=integral(Fp1,-Om_c,0)+integral(Fp1,0,Om_c)+integral(Fm1,-Om_c,0)+integral(Fm1,0,Om_c);
    end
end
%figure(1)
%plot(Omega,Ip,'.-r',Omega,Im,'.-b')
figure(2)
plot(t/Gap,C.^-1,'.-')
hold on;
end


function f=Kubo_0T(G1)
im=sqrt(-1);eta=0.01;t1=100;tau=1/(2*eta);
%following way of defining the paramters are from FW paper but scaled down
Gap=2;%in meV
W=27*Gap;%in meV
m1=1;m2=10*m1;mplus=m2*m1/(m2+m1);mminus=m2*m1/(m2-m1);%in m_e
mu0=W*m2/(m1+m2);%in meV
kconst=6.06*10^-3;
kf=kconst*sqrt(2*mplus*W);%in units of a0^-1
v=Gap/(2*kf);%in units of meV.a0
G2=0.1*Gap;
Gm=G1+G2;gm=G1-G2;

Econst=13605;%hbar^2/(2.me.a0^2)
ek1=@(k) Econst*k.^2/(m1)-mu0;ek2=@(k) Econst*k.^2/(m2)-W+mu0;
delta=@(k) v*(k);
Ek1=@(k) 0.5*(ek1(k)-ek2(k)-im*Gm);
Ek2=@(k) 0.5*sqrt((ek1(k)+ek2(k)-im*gm).^2+4*abs(delta(k).^2));

Ekp=@(k) real(Ek1(k)+Ek2(k));taup=@(k) -1./(2*imag(Ek1(k)+Ek2(k)));
Ekm=@(k) real(Ek1(k)-Ek2(k));taum=@(k) -1./(2*imag(Ek1(k)-Ek2(k)));

vx1=@(k) 0.5*Econst*(2*k/m1-2*k/m2);
vx2=@(k) real(0.25*(2*(ek1(k)+ek2(k)-im*gm).*(Econst*2*k/m1+Econst*2*k/m2)+8*v^2*k)./sqrt((ek1(k)+ek2(k)-im*gm).^2+8*abs(delta(k)).^2));
vxp=@(k) vx1(k)+vx2(k);vxm=@(k) vx1(k)-vx2(k);
dosp=@(k,omega) -1/pi*imag(1./(omega-Ek1(k)-Ek2(k)));
dosm=@(k,omega)-1/pi*imag(1./(omega-Ek1(k)+Ek2(k)));

tht_cont=2/3;

Fp=@(k) tht_cont*k.*vxp(k).^2.*dosp(k,0).^2;
Fm=@(k) tht_cont*k.*vxm(k).^2.*dosm(k,0).^2;
f=integral(Fp,0.0,kf)+integral(Fp,kf,1)+integral(Fm,0.0,kf)+integral(Fm,kf,1);
end
