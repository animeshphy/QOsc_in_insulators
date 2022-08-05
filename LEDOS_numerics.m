%This code does the following
%for different temperatures
%for different combination of gammas
%we calculate dos(T,gm,E,B)
%Do the temperature average to get ledos(T,gm,B) and store it
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
%V2=0;
%definining Matrix G
Ngm=1;
G=zeros(Ngm,2);
G(:,2)=0.0005*Gap*ones(Ngm,1);
G(:,1)=linspace(0.0005,0.0005,Ngm)'*Gap;
%G(Ngm,1)=0.002*Gap;
%G(1,2)=0.002*Gap;

NBinv=1501;b1=(2.5*const)/(m1*mu0);b2=(17.5*const)/(m1*mu0);

Binv=linspace(b1,b2,NBinv);%the range selects the range of n going from 9.5 to 17.5

%temperature in meV
nT=4;
T1=logspace(-2,log(0.4*Gap)/log(10),nT);
T1=[0.035 0.1 0.2 0.4]*Gap;
beta1=T1.^-1;

Nbin=15001;%this is the energy binning
A=cell(1,nT);%DOS_E=A;
for i=1:nT
A{1,i}=zeros(Ngm,NBinv);%this is the DOS as a function of E
%DOS_E{1,i}=zeros(1,Nbin);
end

parfor np=1:nT
    T=T1(np);
    E=mu0+linspace(-10*T,10*T,Nbin);
    for k=1:Ngm
        G1=G(k,1);G2=G(k,2);
        for j=1:NBinv
            DOS_E=zeros(1,Nbin);
            for i=1:Nbin
                n_c=floor(E(i)/(const/(m1*max(Binv)))); %gives the LL index coming from lighter mass
                n_c=2*n_c;
                if n_c<1
                    n_c=1;
                end
                DOS_E(1,i)=DOS(E(i),n_c,Binv(j)^-1,m1,m2,const,G1,G2,V2,W);
            end
            %the following is where we do the temperature averaging
            f=@(x)interp1(E,DOS_E,x);
            beta=T^-1;
            df=@(x) 0.5*beta*exp(-beta*(x-mu0))./(1+exp(-beta*(x-mu0))).^2;
            F=@(x) f(x).*df(x);
            A{1,np}(k,j)=integral(F,min(E),max(E));
        end
    end
end

%Data storing sectio
%save(['ledos_m2_m1_10_1_mult_G1_final_3.mat'],'A','Binv','T1','G')
save(['test_25072022.mat'],'A','Binv','T1','G')
%tocBytes(gcp);
toc
