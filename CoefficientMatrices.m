%function [PhiL,PhiR,vi]=CoefficientMatrices(A,Bs,Ds,taus,cV,p,rmax,d,dt,order,par)
function [PhiL,PhiR,vs]=CoefficientMatrices(systemfun)
A=systemfun.A;
Bs=systemfun.Bs;
Ds=systemfun.Ds;
taus=systemfun.taus;
cV=systemfun.cV;
p=systemfun.p;
rmax=systemfun.rmax;
d=systemfun.d;
dt=systemfun.dt;
par=systemfun.par;


dsize=-(d-1):0;

Ni=length(Bs);
if (all( diff([length(Bs),length(Ds),length(Ds)])))
    error('Dimension missmatch in delay terms')
end

Phi = sparse(1:p*d,1:p*d,1.0,p*d,(p+rmax+1)*d);
vs = zeros((rmax+1)*d,1);
Phi = full(Phi);
E=eye(d,d);

for ti=0:p-1
    t=ti*dt+dt/2;
    P=expm(A(t,par)*dt);
    Phi((p-ti)*d+dsize,(p-ti+1)*d+dsize) = -P;
    for ni=1:Ni
        taui=taus{ni}(t,par);
        rloc=floor(taui/dt);
        R=Bs{ni}(t,par)*Ds{ni}(t/2,par)*dt;
        %R=(P-eye(size(P)))*(A(t,par)\Bs{ni}(t,par)*Ds{ni}(t/2,par))*dt;
        %%R=(P-eye(size(P)))*inv(A(t,par))*Bs{ni}(t,par)*Ds{ni}(t/2,par)*dt;

        Phi((p-ti)*d+dsize,(p-ti+1+rloc)*d+dsize) =-R;

        %Phi((p-ti)*d+dsize,(p-ti+1+rloc)*d+dsize) =-R/2;
        %Phi((p-ti)*d+dsize,(p-ti+1+rloc-1)*d+dsize) =-R/2;
    end
    vs((p-ti)*d+dsize,1)=cV(t,par)*dt;
end


PhiL=Phi(1:p*d,1:p*d);
PhiR=-Phi(1:p*d,p*d+1:end);


Next=(rmax+1-p)*d;%extendson size, size of the overlap
PhiL(end+1:end+Next,end+1:end+Next)=eye(Next,Next);
PhiR(end+1:end+Next,1:Next)=eye(Next,Next);

%  figure(1234)
%  subplot(1,3,1)
%  spy(Phi),grid on
%  subplot(1,3,2)
%  spy(PhiL),grid on
%  subplot(1,3,3)
%  spy(PhiR),grid on