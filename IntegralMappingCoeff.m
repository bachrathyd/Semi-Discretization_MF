function v=IntegralMappingCoeff(s,systemfun)


%M=systemfun.M;
%rlocs=systemfun.rlocs;
%vs=systemfun.vs;

smat=reshape(s,systemfun.d,[]);
smatconti=[smat,zeros(size(smat,1),systemfun.p)];


if systemfun.p~=size(systemfun.M,4)
    error("Size difference in p and M")
end
for ti=0:systemfun.p-1
    tishift=ti+(size(smat,2));%for the proper indexing of smatconti
    %tishift=ti+(rmax);%for the proper indexing of smatconti

    Incement=zeros(systemfun.d,1);
   
    for mi=1:size(systemfun.M,3)
        Incement=Incement+systemfun.M(:,:,mi,ti+1)*...
            smatconti(:,tishift-systemfun.rlocs(mi,ti+1));
    end

    
    smatconti(:,tishift+1) =Incement+ systemfun.vs(:,ti+1);
%     plot(smatconti')
end

v=reshape(smatconti(:,(end-systemfun.rmax):end),[],1);%ez itt inkább R+1