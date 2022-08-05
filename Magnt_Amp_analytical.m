%Following calculation is for finite T magnetization for LF's model
%Here we directly calculate the 
%to calculate DOS for LF parametres
clc;clear all;
tic
im=sqrt(-1);
g_const=1;%2.418e-8;% this is the degeneracy od LLs for 1mmx1mm area
const=1.157*10^-1;%this is e\hbar/m_e  in units of meV*m_e/T
const1=4.25*10^-6;%a0^2*e/(hbar)
Gap=2;%in meV
%mu=25*Gap;%in meV
m1=1;m2=10*m1;mplus=m2*m1/(m2+m1);%mminus=m2*m1/(m2-m1);%in units of mass of electron
W=27*Gap;%mu*2*mminus/mplus;%in meV
mu0=W*m2/(m1+m2);
Ev=Gap^2/(8*W*(mplus/m1));%this is the low energy scale in the problem m1V^2/hbar^2
Eb=5*sqrt(mu0*Ev);%this defines the magnetic field at which we're calculating
%sqrt(Eb*Ev) is the gap scale in this problem=eBhbar/m1
V2=Ev*const/m1;%in meV^2/T
Eb2=Eb*m1/m2;
Ehyb=Eb*Ev;

Ngm=5;
Gm1=linspace(0.0,0.1,Ngm)*Gap;
G2=0.1*Gap;

Ndelmu=7;
Delmu=linspace(0,Eb/2,Ndelmu);%Delmu(Ndelmu)=2*gm;
Tn=100;
T1=logspace(-2,.3,Tn);

%%
% %Following calculation is done using the both symmetric and antisym combination of
% %energies without leaving higher order correction of v
% %For the first solution
M=zeros(Ndelmu,Tn);

for i=1:Ngm
    if i==1
        G1=0;
        G2=0;
    else
    G1=Gm1(i);
    G2=0.1*Gap;
    end
    mu=mu0;
    for nT=1:Tn
    T=T1(nT);
    Chi=2*pi^2*T;%/Eb2;
    ll1=LL2s_full(0,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
    c=Chi*exp(-2*pi*ll1);
    M(i,nT)=c;
    n=0;
    while c>M(i,nT)/2000
        n=n+1
        ll1=LL2s_full(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
        c=Chi*exp(-2*pi*ll1);
        M(i,nT)=M(i,nT)+c;
        if n>300
            break
        end
    end
    end
end
M1=zeros(Ngm,Tn);

for i=1:Ngm
    if i==1
        G1=0;
        G2=0;
    else
    G1=Gm1(i);
    G2=0.1*Gap;
    end
    mu=mu0;
    for nT=1:Tn
    T=T1(nT);
    Chi=2*pi^2*T;%/Eb;
    ll1=LL1s_full(0,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
    c=Chi*exp(-2*pi*ll1);
    M1(i,nT)=c;
    n=0;
    while c>M1(i,nT)/2000
        n=n+1
        ll1=LL1s_full(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
        c=Chi*exp(-2*pi*ll1);
        M1(i,nT)=M1(i,nT)+c;
        if n>300
            break
        end
    end
    end
end

for i=1:Ngm
    if i==1
        G1=0;
        G2=0;
    else
    G1=Gm1(i);
    G2=0.1*Gap;
    end
    mu=mu0;
    for nT=1:Tn
    T=T1(nT);
    Chi=2*pi^2*T;%/Eb2;
    ll1=LL2a_full(0,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
    c=Chi*exp(-2*pi*ll1);
    M(i,nT)=M(i,nT)+c;
    n=0;
    while c>M(i,nT)/2000
        n=n+1
        ll1=LL2a_full(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
        c=Chi*exp(-2*pi*ll1);
        M(i,nT)=M(i,nT)+c;
        if n>300
            break
        end
    end
    end
end

for i=1:Ngm
    if i==1
        G1=0;
        G2=0;
    else
    G1=Gm1(i);
    G2=0.1*Gap;
    end
    mu=mu0;
    for nT=1:Tn
    T=T1(nT);
    Chi=2*pi^2*T;%/Eb;
    ll1=LL1a_full(0,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
    c=Chi*exp(-2*pi*ll1);
    M1(i,nT)=M1(i,nT)+c;
    n=0;
    while c>M1(i,nT)/2000
        n=n+1
        ll1=LL1a_full(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
        c=Chi*exp(-2*pi*ll1);
        M1(i,nT)=M1(i,nT)+c;
        if n>300
            break
        end
    end
    end
end
%M1=2*M1;
% %Plot section
figure(2)
for i=1:Ngm
plot(T1/Gap,(M(i,:)+M1(i,:))/(Eb+Eb2),'.-')
hold on;
end


%%Use the following when you need MvsT as \Gamma changes for E^s under v^2
%%approximation
% %For the first solution
% M=zeros(Ngm,Tn);
% 
% for i=1:Ngm
%     G1=Gm1(i);
%     delmu=Delmu(2);
%     mu=mu0+delmu;
%     for nT=1:Tn
%     T=T1(nT);
%     Chi=2*pi^2*T/Eb2;
%     ll1=LL2(0,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
%     c=Chi*exp(-2*pi*ll1);
%     M(i,nT)=c;
%     n=0;
%     while c>M(i,nT)/2000
%         n=n+1
%         ll1=LL2(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
%         c=Chi*exp(-2*pi*ll1);
%         M(i,nT)=M(i,nT)+c;
%         if n>300
%             break
%         end
%     end
%     end
% end
% M=2*M;
% 
% 
% %for the second solution
% M1=zeros(Ndelmu,Tn);
% 
% for i=1:Ngm
%     G1=Gm1(i);
%     delmu=Delmu(2);
%     mu=mu0+delmu;
%     for nT=1:Tn
%     T=T1(nT);
%     Chi=2*pi^2*T/Eb;
%     ll1=LL1(0,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
%     c=Chi*exp(-2*pi*ll1);
%     M1(i,nT)=c;
%     n=0;
%     while c>M1(i,nT)/2000
%         n=n+1
%         ll1=LL1(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2);
%         c=Chi*exp(-2*pi*ll1);
%         M1(i,nT)=M1(i,nT)+c;
%         if n>300
%             break
%         end
%     end
%     end
% end
% M1=2*M1;
% %Plot section
% figure(2)
% for i=1:Ngm
% plot(T1*10,M(i,:)+M1(i,:),'.-')
% hold on;
% end


%We are interested in the imaginary part of l and the following function
%does that
function a=LL1(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2)
i=sqrt(-1);
omegan=(2*n+1)*pi*T;
a=0;
a=a+i*(omegan+G1)/Eb;
o1=-2*Ehyb*(Eb*(W-i*G2)+i*Eb2*G1-(mu+i*omegan)*(Eb-Eb2))/((mu+i*omegan)*(Eb+Eb2)+i*G1*Eb2+i*G2*Eb-(W+Eb2)*Eb)^2;
a=a+o1*((mu+i*omegan)*(Eb+Eb2)+i*G1*Eb2+i*G2*Eb-(W+Eb2)*Eb)/(2*Eb*Eb2);
a=abs(imag(a));
end

function a=LL2(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2)
i=sqrt(-1);
omegan=(2*n+1)*pi*T;
a=0;
a=a-i*(omegan+G2)/Eb2;
o1=-2*Ehyb*(Eb*(W-i*G2)+i*Eb2*G1-(mu+i*omegan)*(Eb-Eb2))/((mu+i*omegan)*(Eb+Eb2)+i*G1*Eb2+i*G2*Eb-(W+Eb2)*Eb)^2;
a=a-o1*((mu+i*omegan)*(Eb+Eb2)+i*G1*Eb2+i*G2*Eb-(W+Eb2)*Eb)/(2*Eb*Eb2);
a=abs(imag(a));
end
function f=LL1s_full(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2)
i=sqrt(-1);
omegan=(2*n+1)*pi*T;
a=-Eb*Eb2;
b=Eb*(Eb2+W-i*G2)+Eb2*i*G1-(mu+i*omegan)*(Eb-Eb2)-2*Ehyb;
c=-(mu+i*omegan)*(Eb2+W-i*G1-i*G2)-i*G1*(W+Eb2-i*G2)+(mu+i*omegan)^2+Ehyb;
f=-b/(2*a);
f=f+sqrt(b^2-4*a*c)/(2*a);
 f=abs(imag(f));
end
function f=LL2s_full(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2)
i=sqrt(-1);
omegan=(2*n+1)*pi*T;
a=-Eb*Eb2;
b=Eb*(Eb2+W-i*G2)+Eb2*i*G1-(mu+i*omegan)*(Eb-Eb2)-2*Ehyb;
c=-(mu+i*omegan)*(Eb2+W-i*G1-i*G2)-i*G1*(W+Eb2-i*G2)+(mu+i*omegan)^2+Ehyb;
f=-b/(2*a);
f=f-sqrt(b^2-4*a*c)/(2*a);
 f=abs(imag(f));
end
function f=LL1a_full(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2)
i=sqrt(-1);
omegan=(2*n+1)*pi*T;
a=-Eb*Eb2;
b=Eb*(Eb2+W-i*G2)+Eb2*i*G1-(mu+i*omegan)*(Eb-Eb2)-2*Ehyb;
c=-(mu+i*omegan)*(W-Eb-i*G1-i*G2)-(i*G1+Eb)*(W-i*G2)+(mu+i*omegan)^2+Ehyb;
f=-b/(2*a);
f=f-sqrt(b^2-4*a*c)/(2*a);
 f=abs(imag(f));
end
function f=LL2a_full(n,T,Eb,Eb2,Ehyb,mu,W,G1,G2)
i=sqrt(-1);
omegan=(2*n+1)*pi*T;
a=-Eb*Eb2;
b=Eb*(Eb2+W-i*G2)+Eb2*i*G1-(mu+i*omegan)*(Eb-Eb2)-2*Ehyb;
c=-(mu+i*omegan)*(W-Eb-i*G1-i*G2)-(i*G1+Eb)*(W-i*G2)+(mu+i*omegan)^2+Ehyb;
f=-b/(2*a);
f=f+sqrt(b^2-4*a*c)/(2*a);
 f=abs(imag(f));
end