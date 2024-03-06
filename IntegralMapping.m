function v=IntegralMapping(s,systemfun)
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

smat=reshape(s,d,[]);
smatconti=[smat,zeros(size(smat,1),p)];


Ni=length(Bs);
if (all( diff([length(Bs),length(Ds),length(Ds)])))
    error('Dimension missmatch in delay terms')
end


for ti=0:p-1
    t=ti*dt+dt/2;
    tishift=ti+(size(smat,2));%for the proper indexing of smatconti
    %tishift=ti+(rmax);%for the proper indexing of smatconti
    
    
    %P=E+A(t,par)*dt+A(t,par)*A(t,par)*dt^2/2+A(t,par)*A(t,par)*A(t,par)*dt^2/2/3;
    P=expm(A(t,par)*dt);
    %Phi((p-ti)*d+dsize,(p-ti+1)*d+dsize) = -P;
    Incement=P*smatconti(:,tishift);
    

    for ni=1:Ni
        taui=taus{ni}(t,par);
        rloc=floor(taui/dt);
        R=Bs{ni}(t,par)*Ds{ni}(t/2,par)*dt;

        %Phi((p-ti)*d+dsize,(p-ti+1+rloc)*d+dsize) =-R;
        Incement=Incement+R*smatconti(:,tishift-rloc);
    end

    
    %vs((p-ti)*d+dsize,1)=cV(t,par)*dt;
    smatconti(:,tishift+1) =Incement+ cV(t,par)*dt;
%     plot(smatconti')
end

%v=reshape(smatconti(:,(end-p):end),[],1);
v=reshape(smatconti(:,(end-systemfun.rmax):end),[],1);%ez itt inkább R+1