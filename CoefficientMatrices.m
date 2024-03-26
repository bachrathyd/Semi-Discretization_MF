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

m=p*d;
n=(p+rmax+1)*d;

i=(1:p*d)';
j=(1:p*d)';
v=ones(p*d,1);

vs = zeros((rmax+1)*d,1);

E=eye(d,d);

for ti=0:p-1 %all time step 
    t=ti*dt+dt/2;

    %P=expm(A(t,par)*dt); - semi-disc method
    P=E+A(t,par)*dt+A(t,par)*A(t,par)*dt^2/2; %- second order apporximation

    %Phi((p-ti)*d+dsize,(p-ti+1)*d+dsize) = -P;  % filling the full matrixes
    %[I,J]=ndgrid((p-ti)*d+dsize,(p-ti+1)*d+dsize); % generating the correspoind indices %slow version
    I=repmat(((p-ti)*d+dsize)',1,2);% generating the correspoind indices %fast version
    J=repmat((p-ti+1)*d+dsize,2,1);% generating the correspoind indices %fast version
    i=[i;I(:)];
    j=[j;J(:)];
    v=[v;-P(:)];


    for ni=1:Ni
        taui=taus{ni}(t,par);
        rloc=floor(taui/dt);

        %R=(P-eye(size(P)))*(A(t,par)\Bs{ni}(t,par)*Ds{ni}(t/2,par))*dt; %-semi-disc method 0th \\TODO: check
        %%R=(P-eye(size(P)))*inv(A(t,par))*Bs{ni}(t,par)*Ds{ni}(t/2,par)*dt;%- semi-disc method 1st \\TODO: check
        R=Bs{ni}(t,par)*Ds{ni}(t,par)*dt; % 0th Full discretization

        %Phi((p-ti)*d+dsize,(p-ti+1+rloc)*d+dsize) =-R;  % filling the full matrixes
        %[I,J]=ndgrid((p-ti)*d+dsize,(p-ti+1+rloc)*d+dsize); % generating the correspoind indices %slow version

        I=repmat(((p-ti)*d+dsize)',1,2);% generating the correspoind indices %fast version
        J=repmat((p-ti+1+rloc)*d+dsize,2,1);% generating the correspoind indices %fast version
        i=[i;I(:)];
        j=[j;J(:)];
        v=[v;-R(:)];


    end
    vs((p-ti)*d+dsize,1)=cV(t,par)*dt;
end
Phi = sparse(i,j,v,m,n);


PhiL=Phi(1:p*d,1:p*d);
PhiR=-Phi(1:p*d,p*d+1:end);

Next=(rmax+1-p)*d;%extendson size, size of the overlap
PhiL(end+1:end+Next,end+1:end+Next)=eye(Next,Next);
PhiR(end+1:end+Next,1:Next)=eye(Next,Next);

%
% PhiL=sparse(PhiL);
% PhiR=sparse(PhiR);
%  figure(1234)
%  subplot(1,3,1)
%  spy(Phi),grid on
%  subplot(1,3,2)
%  spy(PhiL),grid on
%  subplot(1,3,3)
%  spy(PhiR),grid on