 clc;clear all;

 % %Here I calculate the DOS(0) as a function of magnetic field
% %All the calculation are done in SI units
% %B=Newton.sec/coulomb.metre
% im=sqrt(-1);
% g_const=2.418e-8;% this is the degeneracy od LLs for 1mmx1mm area
% const=1.16*10^-1;%this is e\hbar/m_e  in units of meV*m_e/T
% const1=4.25*10^-6;%a0^2*e/(hbar)
% hop=640;%in meV
% m1=13605/hop;m2=50*m1;
% mplus=m2*m1/(m2+m1);mminus=m2*m1/(m2-m1);
% W=.5*hop;
% mu0=W*m2/(m1+m2);%in meV
% V2=(0.015*hop)^2*const1;%in meV^2
% Ev=V2*m1/const;%this is the low energy scale in the problem mv^2/hbar^2
% Gap=sqrt(8*W*Ev);

%LF parametres
im=sqrt(-1);
g_const=2.418e-8;% this is the degeneracy od LLs for 1mmx1mm area
const=1.16*10^-1;%this is e\hbar/m_e  in units of meV*m_e/T
const1=4.25*10^-6;%a0^2*e/(hbar)
Gap=2;
%mu=25*Gap;
m1=1;m2=10*m1;
Gap1=Gap*sqrt(4*m1*m2)/(m1+m2);
mplus=m2*m1/(m2+m1);%mminus=m2*m1/(m2-m1);
W=27*Gap;%mu*2*mminus/mplus;
mu0=W*m2/(m1+m2);%in meV
% mu0=50*Gap;%mu*2*mminus/mplus;
% W=mu0*(m1+m2)/m2;%in meV
Ev=Gap^2/(8*W*(mplus/m1));%this is the low energy scale of the problem mv^2/hbar^2
V2=Ev*const/m1;
perd2=0;

sz=1000;sz1=30;
T1=logspace(-2,2,sz1);%temp is in meV
T=1;%in meV
Binverse=linspace(0.4,25,sz)*const/(mu0*m1);E=zeros(1,sz);
k=1;
%In the following section we plot the LEDOS as a function of Binv
x1=cell(1,70);y1=x1;y2=x1;x2=x1;x3=x1;y3=x1;x4=x1;y4=x1;
for n=1:1:70
    count=ones(1,4);
for i=1:sz
    b=Binverse(i)^-1;
    %nf=floor((-mu-g/2)*m1/(b*const));%this correspods to the Landau level around the fermi energy
        en1= (n-0.5)*const*b/m1;e1=(n+0.5)*const*b/m1;
        en2= W-(n+0.5)*const*b/m2;e2=W-(n-0.5)*const*b/m2;
        Enp=0.5*(en1+en2+sqrt((en1-en2)^2+8*n*V2*b))-mu0;
        Enp1=0.5*(e1+e2+sqrt((e1-e2)^2+8*n*V2*b))-mu0;
        Enm=0.5*(en1+en2-sqrt((en1-en2)^2+8*n*V2*b))-mu0;
        Enm1=0.5*(e1+e2-sqrt((e1-e2)^2+8*n*V2*b))-mu0;
        E_cut=W/8;
        
        if abs(Enp1)<E_cut
        y1{1,n}(count(1))=Enp1;
        x1{1,n}(count(1))=b^-1;
        count(1)=count(1)+1;
        end
        if abs(Enm1)<E_cut
        y2{1,n}(count(2))=Enm1;
        x2{1,n}(count(2))=b^-1;
        count(2)=count(2)+1;
        end
        if abs(Enp)<E_cut
        y3{1,n}(count(3))=Enp;
        x3{1,n}(count(3))=b^-1;
        count(3)=count(3)+1;
        end
        
        if abs(Enm)<E_cut
        y4{1,n}(count(4))=Enm;
        x4{1,n}(count(4))=b^-1;
        count(4)=count(4)+1;
        end

end
%[a,b]=min(y3{1,n});
 %A(1,k)=a;A(2,k)=Binverse(b)^-1;
 k=k+1;
end
c1=[0,0.447,0.741];
c2=[0.85,0.325,0.098];
c3=[0.929,0.694,0.125];
c4=[0.494,0.184,0.556]

for n=1:70
plot(m1*x3{1,n}*mu0/const,y3{1,n},'.-','color',c1)
hold on
end
for n=1:70
    plot(m1*x1{1,n}*mu0/const,y1{1,n},'.-','color',c2)
hold on
end

for n=1:70
plot(m1*x2{1,n}*mu0/const,y2{1,n},'.-','color',c3)
 hold on
end
for n=1:70
%if count(4) ~= 1
plot(m1*x4{1,n}*mu0/const,y4{1,n},'.-','color',c4)
hold on
%end
end
%

%In the following section we plot the periodicities coming from E_+ and E_-
Nx=21;
x=linspace(1,Nx,Nx);
y=x;
for n=1:21
    [a,b]=max(y2{1,n});
    a0=x2{1,n}(b);
    [a,b]=max(y2{1,n+1});
    b0=x2{1,n+1}(b);
    y(n)=(b0-a0)*m1*mu0/const;
end
figure(5)
plot(x,y,'.-r')

Nx=21;
x=linspace(1,Nx,Nx);
y=x;
for n=1:21
    [a,b]=min(y3{1,n});
    a0=x3{1,n}(b);
    [a,b]=min(y3{1,n+1});
    b0=x3{1,n+1}(b);
    y(n)=(b0-a0)*m1*mu0/const;
end
figure(6)
plot(x,y,'.-b')

