%This is a test code for the main temperature averaging calculation
clc;clear all;
tic
%Initialisation
im=sqrt(-1);
g_const=1;%2.418e-8;% this is the degeneracy od LLs for 1mmx1mm area
const=1.157*10^-1;%this is e\hbar/m_e  in units of meV*m_e/T
const1=4.25*10^-6;%a0^2*e/(hbar)
Gap=2;%in meV
m1=1;m2=10*m1;mplus=m2*m1/(m2+m1);mminus=m2*m1/(m2-m1);%in units of mass of electron
W=27*Gap;%in meV
mu0=W*m2/(m1+m2);
Ev=Gap^2/(8*W*(mplus/m1));%this is the low energy scale in the problem m1V^2/hbar^2
V2=Ev*const/m1;%in meV^2/T
%Ev=0;
%definining Matrix G
Ngm=6;
G=zeros(Ngm,2);
G(1:Ngm,2)=0.05*Gap*ones(Ngm,1);
G(1:Ngm,1)=linspace(0.05,0.55,Ngm)'*Gap;
%G(1,1)=0.0005*Gap;
%G(1,2)=0.0005*Gap;

NBinv=400;b1=(9.5*const)/(m1*mu0);b2=(17.5*const)/(m1*mu0);

Binv=linspace(b1,b2,NBinv);%the range selects the range of n going from 9.5 to 17.5

%temperature in meV
nT=60;
T1=logspace(-2,log(2*Gap)/log(10),nT);
beta1=T1.^-1;

Nbin=4001;%this is the energy binning
A=cell(1,nT);%DOS_E=A;
for i=1:nT
A{1,i}=zeros(Ngm,NBinv);%this is the DOS as a function of E
%DOS_E{1,i}=zeros(1,Nbin);
end
A1=A;A2=A;
parfor np=1:nT
    T=T1(np);
    E=mu0+linspace(-10*T,10*T,Nbin);
    for k=1:Ngm
        G1=G(k,1);G2=G(k,2);
        for j=1:NBinv
            DOS_E=zeros(2,Nbin);
            for i=1:Nbin
                [f1,f2]=DOS_LF_full_1(E(i),Binv(j)^-1,m1,m2,const,G1,G2,W,Gap,Ev);
                DOS_E(1,i)=f1;
                DOS_E(2,i)=f2;
            end
            %the following is where we do the temperature averaging
            beta=T^-1;
            df=@(x) 0.5*beta*exp(-beta*(x-mu0))./(1+exp(-beta*(x-mu0))).^2;
            f=@(x)interp1(E,DOS_E(1,:),x);
            F=@(x) f(x).*df(x);
            A1{1,np}(k,j)=integral(F,min(E),max(E));
            f=@(x)interp1(E,DOS_E(2,:),x);
            F=@(x) f(x).*df(x);
            A2{1,np}(k,j)=integral(F,min(E),max(E));
            A{1,np}(k,j)=A1{1,np}(k,j)+A2{1,np}(k,j);
        end
    end
end

%Data storing sectio
save(['ledos_analytics_final_m1_m2_10_1_mult_G1_G2_pt05.mat'],'A','A1','A2','Binv','T1','G')
%tocBytes(gcp);
toc