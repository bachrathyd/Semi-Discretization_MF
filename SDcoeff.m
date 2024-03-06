function systemfun=SDcoeff(systemfun)
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



Ni=length(Bs);
if (all( diff([length(Bs),length(Ds),length(Ds)])))
    error('Dimension missmatch in delay terms')
end

M=zeros(d,d,Ni+1,p);
rlocs=zeros(Ni+1,p);
vs=zeros(d,p);


for ti=0:p-1
    t=ti*dt+dt/2;
    tishift=ti+1;%for the proper indexing of smatconti
    
    P=expm(A(t,par)*dt);
    M(:,:,1,tishift)=P;
    rlocs(1,tishift)=0;
    

    for ni=1:Ni
        taui=taus{ni}(t,par);
        rloc=floor(taui/dt);
        R=Bs{ni}(t,par)*Ds{ni}(t/2,par)*dt;
        M(:,:,1+ni,tishift)=R;
        rlocs(1+ni,tishift)=rloc;
    end
    
    vs(:,tishift)=cV(t,par)*dt;

end

systemfun.M=M;
systemfun.rlocs=rlocs;
systemfun.vs=vs;