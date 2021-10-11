function M1=LogCinit(sinterdiff,T,W,isclose)
% *** Corrections to Logarithmic potential log(|tau-z|) ***   
% *** block-tri-diagonal output ***
% *** iper=0,1 (0 is open arc, 1 is closed contour) ***
npan = length(sinterdiff);
ngl = length(T);
N = ngl*npan;
M1=zeros(N);
A=ones(ngl);
for k=2:ngl
  A(:,k)=T.*A(:,k-1);
end
if isclose
  kstart=1;
  kend=npan;  
else
  kstart=2;
  kend=npan-1;  
end
% *** central blocks ***
TMP=-log(abs(T-T'));
TMP(1:ngl+1:ngl^2)=0;
TMP=TMP+chnk.quadjh.WfrakLinit(A,0,1,T)./W(:,ones(1,ngl))';
for k=1:npan
  myind=(k-1)*ngl+1:k*ngl;
  M1(myind,myind)=TMP;
end
% *** superdiagonal blocks (targets to the left) ***
for k=kstart:npan
  myinds=(k-1)*ngl+1:k*ngl;
  km1=mod(k-2,npan)+1;
  mi=(km1-1)*ngl+1:km1*ngl;       
  alpha=sinterdiff(km1)/sinterdiff(k);          
  [TMP,accept,na]=chnk.quadjh.WfrakLinit(A,-1-alpha,alpha,T);
  mi=mi(accept);
  for nj=1:ngl
    M1(mi,myinds(nj))=-log(abs(T(nj)+1+(1-T(accept))*alpha));       
  end
  M1(mi,myinds)=M1(mi,myinds)+TMP./W(:,ones(1,na))';
end
% *** subdiagonal blocks (targets to the right) ***
for k=1:kend
  myinds=(k-1)*ngl+1:k*ngl;
  kp1=mod(k,npan)+1;
  mi=(kp1-1)*ngl+1:kp1*ngl;       
  alpha=sinterdiff(kp1)/sinterdiff(k);               
  [TMP,accept,na]=chnk.quadjh.WfrakLinit(A,1+alpha,alpha,T);
  mi=mi(accept);
  for nj=1:ngl
    M1(mi,myinds(nj))=-log(abs(T(nj)-1-(T(accept)+1)*alpha));   
  end
  M1(mi,myinds)=M1(mi,myinds)+TMP./W(:,ones(1,na))';
end
if isclose
  M1=sparse(M1);
end
