function H=Foo_eig_MDBM(ax,systemfun)

%% version 1: simple calculate the function value(s) for each parameter points
H=zeros(1,size(ax,2));
for k=1:size(ax,2)
    systemfun.par.delta=ax(1,k);
    systemfun.par.b0=ax(2,k);
    if size(ax,1)>2
    systemfun.par.eps=ax(3,k);
    end

    %FASTEST %e.g: 19s
    [PhiL,PhiR,vi]=CoefficientMatrices(systemfun);
    H(1,k)=log(max(abs(eigs(PhiR,PhiL))));


% %SLIGHLY SLOWER %e.g: 47s
%    systemfun.rmax=ceil(systemfun.par.taumax/systemfun.dt);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);
%     N=(systemfun.rmax+1)*systemfun.d;
%     systemfun=SDcoeff(systemfun);
%     s0=rand(N,1);
%     v0=IntegralMappingCoeff(s0,systemfun);
%     AffineMappingPerturbe=@(s) IntegralMappingCoeff(s+s0,systemfun)-v0;
%     H(1,k)=log(max(abs(eigs(AffineMappingPerturbe,N))));


%     %SLOWEST %e.g: 229s
%     systemfun.rmax=ceil(systemfun.par.taumax/systemfun.dt);% stepsize for the delay %TODO: it can be reduced to r=ceil(par.taumax/dt);
%     N=(systemfun.rmax+1)*systemfun.d;
%     s0=rand(N,1);
%     v0=IntegralMapping(s0,systemfun);
%     AffineMappingPerturbe=@(s) IntegralMapping(s+s0,systemfun)-v0;
%     H(1,k)=log(max(abs(eigs(AffineMappingPerturbe,N))));

end
