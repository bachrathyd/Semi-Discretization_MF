

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
par.T=0.7*pi; %time period of the system %overlap (T > tau)
par.T=2.0*pi; %time period of the system %equal (T = tau)
% par.T=7.0*pi; %time period of the system %missing elements (T > tau)

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

% % % tic
% % % [PhiL,PhiR,vs]=CoefficientMatrices(systemfun);
% % % eigs(PhiR,PhiL)
% % % toc
% % %
% % % tic
% % % eigs(PhiL\PhiR)%ez gyorsabbnak tűnik
% % % toc
% % %
% % % tic
% % % fp=(PhiL-PhiR)\vs;
% % % toc
% % % figure(32)
% % % plot(fp)
% % %
% % % %% Integration - onthe fly cacluation, without storing
% % % % Important if T>tau
% % % rmax=ceil(par.taumax/dt);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);
% % %
% % % systemfun.rmax=rmax;
% % % % rmax=max(ceil(par.taumax/dt),p-1);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);
% % % N=(rmax+1)*d;
% % % v0=rand(N,1);
% % %
% % % % v0LR=PhiL\PhiR*s0
% % %
% % % figure(124),clf
% % % plot((1:(rmax+1)*d),v0), hold on
% % %
% % % for k=1:20
% % %     v0=IntegralMapping(v0,systemfun);
% % %     plot((1:(rmax+1)*d) +k*(p+0.1)*d,v0), hold on % 0.1 is to prevent the overlap
% % %     %plot(v0LR)
% % %     % pause
% % % end
% % %
% % %
% % % figure(32),clf
% % % plot(fp), hold on
% % % plot(v0(end:-1:1)), hold on %the values are stroed in reverse order
% % %
% % % tic
% % % s0=rand(N,1);
% % % v0=IntegralMapping(s0,systemfun);
% % % AffineMappingPerturbe=@(s) IntegralMapping(s+s0,systemfun)-v0;
% % % %LinMappingPerturbe=@(s) IntegralMapping(s,systemfun)
% % % eigs(AffineMappingPerturbe,N)
% % % toc
% % %
% % %
% % % %% Integration based with precomputed coefficient matrices
% % %
% % %
% % % systemfun=SDcoeff(systemfun)
% % % v0C=IntegralMappingCoeff(s0,systemfun);
% % %
% % % figure(54),clf
% % % plot(v0C), hold on
% % % plot(v0,'--',LineWidth=3)
% % %
% % %
% % % tic
% % % AffineMappingPerturbe=@(s) IntegralMappingCoeff(s+s0,systemfun)-v0;
% % % %LinMappingPerturbe=@(s) IntegralMapping(s,systemfun)
% % % eigs(AffineMappingPerturbe,N)
% % % toc
% % %


%% - Time complextiy check -



par=[];
par.delta=3.0;
par.eps=2.0;
par.b0=-0.15;
par.kappa=0.05;
par.T=2*pi; %time period of the system
par.taumax=2*pi; %maximal delay, used for the


% pv=ceil(10.^(1.0:0.1:4.0));
pv=ceil(10.^(1.5:0.2:4.0));
pv=ceil(10.^linspace(1.1,4.5,500));
pv=ceil(10.^linspace(2,5,50));
Tcpu_int=nan*pv;
Tcpu_PhiLR=nan*pv;
Tcpu_Phi=nan*pv;
Neig=15;
mus=nan(Neig,length(pv));
for kp=1:length(pv)
    p=pv(kp)


    dt=par.T/p; %stepsize
    rmax=ceil(par.taumax/dt);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);

    d=size(A(0.0,par),1);%dimension of the sytem (states)

    systemfun.p=p;
    systemfun.rmax=rmax;
    systemfun.d=d;
    systemfun.dt=dt;
    systemfun.par=par;

    N=(rmax+1)*d;
    s0=rand(N,1);

    % ~linear in time p>~200  %Constant multiplier depends on Neig (almost proporionally if Neig larger (>10) )

    systemfun=SDcoeff(systemfun);
    v0=IntegralMappingCoeff(s0,systemfun);
     if p<20000
    tic
    AffineMappingPerturbe=@(s) IntegralMappingCoeff(s+s0,systemfun)-v0;
    %LinMappingPerturbe=@(s) IntegralMapping(s,systemfun)
    mus(:,kp)=eigs(AffineMappingPerturbe,N,Neig);
    Tcpu_int(kp)=toc;%/Neig;
     end
    % ~quadratic in time if p>~200
    rmax=max(ceil(par.taumax/dt),p-1);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);
    systemfun.rmax=rmax;

    tic
    [PhiL,PhiR,vs]=CoefficientMatrices(systemfun);
    mus(:,kp)=eigs(PhiR,PhiL,Neig);
    Tcpu_PhiLR(kp)=toc;
     if p<2000


        % [PhiL,PhiR,vs]=CoefficientMatrices(systemfun);
        tic
        [PhiL,PhiR,vs]=CoefficientMatrices(systemfun);
        mus(:,kp)=eigs(PhiL\PhiR,Neig);
        Tcpu_Phi(kp)=toc;

     end
end

figure(432)
clf
% subplot(2,1,1)
% loglog(pv,abs(mus)),hold on
% subplot(2,1,2)

coefficients = polyfit(log(pv),log(Tcpu_int), 1);
% loglog(pv,Tcpu_int,'DisplayName',num2str(coefficients(1))), hold on
loglog(pv,Tcpu_int,'DisplayName',['Tcpu Int-',num2str(coefficients(1))]), hold on
coefficients = polyfit(log(pv(~isnan(Tcpu_PhiLR))),log(Tcpu_PhiLR(~isnan(Tcpu_PhiLR))), 1);
% loglog(pv,Tcpu_PhiLR,'DisplayName',num2str(coefficients(1))), hold on
loglog(pv,Tcpu_PhiLR,'DisplayName',['Tcpu PhiLR-',num2str(coefficients(1))]), hold on
coefficients = polyfit(log(pv(~isnan(Tcpu_Phi))),log(Tcpu_Phi(~isnan(Tcpu_Phi))), 1);
% loglog(pv,Tcpu_Phi,'DisplayName',num2str(coefficients(1))), hold on
loglog(pv,Tcpu_Phi,'DisplayName',['Tcpu Phi-',num2str(coefficients(1))]), hold on
legend
grid on

axis([10,1e6,1e-5,20])


%% ------------MDBM-----------------
tic
addpath('C:\Users\Bacharthy\Documents\GitHub\MBDM\code_folder')

par=[];
par.delta=nan;
par.eps=1.0;
par.b0=nan;
par.kappa=0.05;
par.T=2*pi; %time period of the system
par.taumax=2*pi; %maximal delay, used for the resolution of the timedelay


p=200; %time steps for a full period
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
ax(1).val=-1:0.5:5;%delta
ax(2).val=-1.5:0.2:1.5;%b0
ax(3).val=linspace(-0.01,5,5);%epsilon
%Foo_eig_MDBM([1.1;3.3],systemfun)
mdbm_options=mdbmset('timelimit',inf);
mdbm_sol=mdbm(ax,'Foo_eig_MDBM',4,mdbm_options,systemfun);
figure(2)
clf
hold on
plotobject=plot_mdbm(mdbm_sol);
set(plotobject,'LineWidth',2)
% view(2)
% xlim([-0.35,5])
% ylim([-1.,1.2])
shading interp
light('Position',[-10 0 0],'Style','local')
light('Position',[0 10 0],'Style','local')
light('Position',[0 0 10],'Style','local')
toc
