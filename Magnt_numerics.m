%Following calculation is for finite T magnetization for FW's model
clc;clear all;

%Fa Wang thing
im=sqrt(-1);
g_const=2.418e-8;% this is the degeneracy od LLs for 1mmx1mm area
const=1.16*10^-1;%this is e\hbar/m_e  in units of meV*m_e/T
const1=4.25*10^-6;%a0^2*e/(hbar)
Gap=2;
m1=1;m2=10*m1;
Gap1=Gap*sqrt(4*m1*m2)/(m1+m2);
mplus=m2*m1/(m2+m1);%mminus=m2*m1/(m2-m1);
W=27*Gap;%mu*2*mminus/mplus;
mu0=W*m2/(m1+m2);%in meV
% mu0=50*Gap;%mu*2*mminus/mplus;
% W=mu0*(m1+m2)/m2;%in meV
Ev=Gap^2/(8*W*(mplus/m1));%this is the low energy scale of the problem mv^2/hbar^2
V2=Ev*const/m1;
V2=0;

sz=1000;
Binv=linspace(0.5,17.5,sz)*const/(mu0*m1);
% %In the following section we plot the QOsc Amplitude as a function of temp
omega=0;
B1=Binv.^-1+.1;
Tn=1;
T1=linspace(.01,Gap,Tn);
Pot=zeros(Tn,sz);Pot1=Pot;
M=Pot;


for i=1:sz
    b=Binv(i)^-1;
    %THIS FOLLOWING CALCULATION DOESN'T NEED TO BE DONE AT EVERY TEMPERATURE
    [Omega_Mu,Mu1]=Omega_T0(b,max(T1),m1,m2,mu0,W,V2,const);
    for nT=1:Tn
    T=T1(nT);
    beta=1/T;
    dfp=@(Mu) beta*exp(-beta*Mu)./(1+exp(-beta*Mu)).^2;
    dfm=@(Mu) beta*exp(beta*Mu)./(1+exp(beta*Mu)).^2;
    Ip=@(Mu) interp1(Mu1,Omega_Mu,Mu).*dfp(Mu);
    Im=@(Mu) interp1(Mu1,Omega_Mu,Mu).*dfm(Mu);
    Muc=10*T;
    Pot(nT,i)=integral(Ip,-Muc,0)+integral(Im,0,Muc);
    end
end

for i=1:sz
    b1=B1(i);
    %THIS FOLLOWING CALCULATION DOESN'T NEED TO BE DONE AT EVERY TEMPERATURE
    [Omega_Mu,Mu1]=Omega_T0(b1,max(T1),m1,m2,mu0,W,V2,const);
    for nT=1:Tn
    T=T1(nT);
    beta=1/T;
    dfp=@(Mu) beta*exp(-beta*Mu)./(1+exp(-beta*Mu)).^2;
    dfm=@(Mu) beta*exp(beta*Mu)./(1+exp(beta*Mu)).^2;
    Ip=@(Mu) interp1(Mu1,Omega_Mu,Mu).*dfp(Mu);
    Im=@(Mu) interp1(Mu1,Omega_Mu,Mu).*dfm(Mu);
    Muc=10*T;
    Pot1(nT,i)=integral(Ip,-Muc,0)+integral(Im,0,Muc);
    end
end

%Plot section
M=Pot1-Pot;
for i=1:Tn
M(i,:)=M(i,:)./(Binv.^-1-B1);
end
figure(1)
for i=1:Tn
plot(m1*Binv*mu0/const,M(i,:),'.-')
hold on;
end
figure(2)
for i=1:Tn
plot(m1*Binv*mu0/const,Pot(i,:),'.-')
hold on;
end

%Following function defines the energy range from given temperature and
%calculates the Grand canonical potential in that range
function [f,Mu]=Omega_T0(B,T,m1,m2,mu,g,V2,const)
sz=3000;%%this is the binning of chemical potential range
f=zeros(1,sz);
Mu=linspace(-10*T,10*T,sz);
nf=floor((-mu-g/2)*m1/(B*const)-1/2);%this correspods to the Landau level around the fermi energy
    nc=500;
    for Mun=1:sz
        Pot1=0;
    for n=0:nc
        if n==0
            E1=ed(n,B,m1,mu,g,const)-Mu(Mun);
            E2=-ef(n,B,m2,mu,g,const)-Mu(Mun);
            Pot1=Pot1-E2*B;
            if E1<0
                Pot1=Pot1+E1*B;
            end
            if E2<0
                Pot1=Pot1+E2*B;
            end
        else
        e1=ed(n,B,m1,mu,g,const);en1=ed(n-1,B,m1,mu,g,const);
        e2=ef(n-1,B,m2,mu,g,const); en2=ef(n,B,m2,mu,g,const);
        E1=0.5*(e1-e2-sqrt((e1+e2)^2+8*n*V2*B))-Mu(Mun);E3=0.5*(en1-en2-sqrt((en1+en2)^2+8*n*V2*B))-Mu(Mun);
        E2=0.5*(e1-e2+sqrt((e1+e2)^2+8*n*V2*B))-Mu(Mun);E4=0.5*(en1-en2+sqrt((en1+en2)^2+8*n*V2*B))-Mu(Mun);
        Pot1=Pot1+e2*B+en2*B;
        if E1<0
            Pot1=Pot1+E1*B;
        end
        if E2<0
            Pot1=Pot1+E2*B;
        end
         if E3<0
            Pot1=Pot1+E3*B;
        end
        if E4<0
            Pot1=Pot1+E4*B;
        end
        end
    end
    f(Mun)=Pot1;
    end
end
function f=ed(n,b,m1,mu0,W,const)
f=(n+0.5)*const*b/m1-mu0;
end
function f=ef(n,b,m2,mu0,W,const)
f=(n+0.5)*const*b/m2-W+mu0;
end