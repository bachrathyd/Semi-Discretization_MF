

%problem definition
A=@(t,par) [0. 1.;...
    -par.delta-par.eps*cos(2*pi/par.T*t) -2*par.kappa];
%B=@(t,par) [0. 0.;...
%    par.b0 0.];
B1=@(t,par) [0.;...
    par.b0];
D1=@(t,par) [1.0, 0.];%the possible improvement is not yet utilized, it is just a preparation for further improvement
tau1=@(t,par) 2*pi ;

cVec = @(t,par) [0.;sin(4*pi/par.T*t)];
%cVec = @(t,par) [0.;0.];

%xp(t)=A(t)*x(t)+B(t)*D(t)*x(t-tau(t)) +c(t)



par=[];
par.delta=3.0;
par.eps=2.0;
par.b0=-0.15;
par.kappa=0.05;
par.T=2*pi; %time period of the system

par.taumax=2*pi; %maximal delay, used for the resolution of the timedelay

p=50; %time steps for a full period
dt=par.T/p; %stepsize
rmax=max(ceil(par.taumax/dt),p-1);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);

d=size(A(0.0,par),1);%dimension of the sytem (states)




systemfun.A=A;
systemfun.Bs={B1};
systemfun.Ds={D1};
systemfun.taus={tau1};
systemfun.cV=cVec;

systemfun.p=p;
systemfun.rmax=rmax;
systemfun.d=d;
systemfun.dt=dt;
systemfun.par=par;

tic
[PhiL,PhiR,vs]=CoefficientMatrices(systemfun);
eigs(PhiR,PhiL)
toc

tic
eigs(PhiL\PhiR)%ez gyorsabbnak tűnik
toc

tic
fp=(PhiL-PhiR)\vs;
toc
figure(32)
plot(fp)

%% Integration - onthe fly cacluation, without storing
N=(rmax+1)*d;
v0=rand(N,1);

% v0LR=PhiL\PhiR*s0

figure(124),clf
plot((1:(rmax+1)*d),v0), hold on

for k=1:20
    v0=IntegralMapping(v0,systemfun);
    plot((1:(rmax+1)*d) +k*(rmax+1)*d,v0), hold on
    %plot(v0LR)
    % pause
end


figure(32),clf
plot(fp), hold on
plot(v0(end:-1:1)), hold on %the values are stroed in reverse order

tic
s0=rand(N,1);
v0=IntegralMapping(s0,systemfun);
AffineMappingPerturbe=@(s) IntegralMapping(s+s0,systemfun)-v0;
%LinMappingPerturbe=@(s) IntegralMapping(s,systemfun)
eigs(AffineMappingPerturbe,N)
toc


%% Integration based with precomputed coefficient matrices


systemfun=SDcoeff(systemfun)
v0C=IntegralMappingCoeff(s0,systemfun);

figure(54),clf
plot(v0C), hold on
plot(v0,'--',LineWidth=3)


tic
AffineMappingPerturbe=@(s) IntegralMappingCoeff(s+s0,systemfun)-v0;
%LinMappingPerturbe=@(s) IntegralMapping(s,systemfun)
eigs(AffineMappingPerturbe,N)
toc



%% - Time complextiy check -



par=[];
par.delta=3.0;
par.eps=2.0;
par.b0=-0.15;
par.kappa=0.05;
par.T=2*pi; %time period of the system
par.taumax=2*pi; %maximal delay, used for the


pv=ceil(10.^(2.3:0.1:4.0));
pv=ceil(10.^(1.0:0.1:3.3));
Tcpu_int=nan*pv;
Tcpu_PhiLR=nan*pv;
Tcpu_Phi=nan*pv;
Neig=1;
mus=nan(Neig,length(pv));
for kp=1:length(pv)
    p=pv(kp)


dt=par.T/p; %stepsize
rmax=max(ceil(par.taumax/dt),p-1);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);

d=size(A(0.0,par),1);%dimension of the sytem (states)

    systemfun.p=p;
    systemfun.rmax=rmax;
    systemfun.d=d;
    systemfun.dt=dt;
    systemfun.par=par;

    N=(rmax+1)*d;
    s0=rand(N,1);
    systemfun=SDcoeff(systemfun);
    v0=IntegralMappingCoeff(s0,systemfun);

       % ~linear in time p>~200
    tic
    AffineMappingPerturbe=@(s) IntegralMappingCoeff(s+s0,systemfun)-v0;
    %LinMappingPerturbe=@(s) IntegralMapping(s,systemfun)
    mus(:,kp)=eigs(AffineMappingPerturbe,N,Neig);
    Tcpu_int(kp)=toc;

    % ~quadratic in time if p>~200
tic
[PhiL,PhiR,vs]=CoefficientMatrices(systemfun);
   mus(:,kp)=eigs(PhiL\PhiR,Neig);
Tcpu_Phi(kp)=toc;
tic
[PhiL,PhiR,vs]=CoefficientMatrices(systemfun);
   mus(:,kp)=eigs(PhiR,PhiL,Neig);
Tcpu_PhiLR(kp)=toc;
end

figure(432)
% % clf
% subplot(2,1,1)
% loglog(pv,abs(mus)),hold on
% subplot(2,1,2)

coefficients = polyfit(log(pv),log(Tcpu_int), 1);
loglog(pv,Tcpu_int,'DisplayName',num2str(coefficients(1))), hold on
coefficients = polyfit(log(pv),log(Tcpu_PhiLR), 1);
loglog(pv,Tcpu_PhiLR,'DisplayName',num2str(coefficients(1))), hold on
coefficients = polyfit(log(pv),log(Tcpu_Phi), 1);
loglog(pv,Tcpu_Phi,'DisplayName',num2str(coefficients(1))), hold on
legend
grid on


%% ------------MDBM-----------------
addpath('C:\Users\Bacharthy\Documents\GitHub\MBDM\code_folder')

par=[];
par.delta=nan;
par.eps=1.0;
par.b0=nan;
par.kappa=0.05;
par.T=2*pi; %time period of the system
par.taumax=2*pi; %maximal delay, used for the resolution of the timedelay


p=130; %time steps for a full period
dt=par.T/p; %stepsize
rmax=max(ceil(par.taumax/dt),p-1);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);



systemfun.A=A;
systemfun.Bs={B1};
systemfun.Ds={D1};
systemfun.taus={tau1};
systemfun.cV=cVec;

systemfun.p=p;
systemfun.rmax=rmax;
systemfun.d=d;
systemfun.dt=dt;
systemfun.par=par;

%the limits and the initial mesh
ax=[];
ax(1).val=-1:0.2:5;%delta
ax(2).val=-1.5:0.2:1.5;%b0
%Foo_eig_MDBM([1.1;3.3],systemfun)
mdbm_options=mdbmset('timelimit',500);
mdbm_sol=mdbm(ax,'Foo_eig_MDBM',2,mdbm_options,systemfun);
figure(2)
% clf
hold on
plotobject=plot_mdbm(mdbm_sol,'k');
set(plotobject,'LineWidth',2)
view(2)
xlim([-0.35,5])
ylim([-1.,1.2])
