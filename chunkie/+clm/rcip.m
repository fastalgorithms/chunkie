  
  function [Pbc,PWbc,starL,circL,starS,circS] = setupRCIP(ngl,ndim,nedge,isstart) 
  % setup for the RCIP forward recursion
  % inputs:
  % ngl - number of Gauss-Legendre nodes
  % nedge - number of edges meeting at the corner
  % isstart - array of size nedge
  %           0 - means the edge ends at the corner w.r.t. the
  %           parameterization of the curve
  %           1 - means the edge starts from the corner w.r.t the
  %           parameterization of the curve
  % ndim - number of equations
  %
  % outputs:
  % Pbc - prolongation matrix
  % PWbc - weighted prolongation matrix
  % starL, circL - bad and good indices for the local system matrix
  % starS, circS - bad and good indices for the preconditioner R
  %
  [T,W] = lege.exps(ngl);  
  [IP,IPW]=IPinit(T,W);
  Pbc = Pbcinit(IP,nedge,ndim);
  PWbc = Pbcinit(IPW,nedge,ndim);
  nedge = 3; % number of edges, nedge=3 means a triple junction
  ndim = 2; % number of equations
  
  % starL - bad indices for the system matrix M
  % circL - good indices for the system matrix M
  starL = [];
  circL = [];
  indg1 = 2*ngl + (1:ngl);
  indb1 = 1:2*ngl;
  
  indg0 = 1:ngl;
  indb0 = ngl + (1:2*ngl);
  
  for i=1:nedge
    if isstart(i) 
      starL = [starL indb1+3*(i-1)*ngl];
      circL = [circL indg1+3*(i-1)*ngl];
    else
      starL = [starL indb0+3*(i-1)*ngl];
      circL = [circL indg0+3*(i-1)*ngl]; 
    end
  end
  
  for i=1:ndim
    starL = [starL starL+3*nedge*ngl]; % type-b mesh has 3 chunks
    circL = [circL circL+3*nedge*ngl];
  end
  
  % starS - bad indices for the preconditioner R
  % circS - good indices for the preconditioner R
  starS = [];
  circS = [];
  indg1 = ngl + (1:ngl);
  indb1 = 1:ngl;
  
  indg0 = 1:ngl;
  indb0 = ngl + (1:ngl);
  
  for i=1:nedge
    if isstart(i)
      starS = [starS indb1+2*(i-1)*ngl]; % type-c mesh has 2 chunks
      circS = [circS indg1+2*(i-1)*ngl];
    else
      starS = [starS indb1+2*(i-1)*ngl];
      circS = [circS indg1+2*(i-1)*ngl]; 
    end
  end
  
  for i=1:ndim
    starS = [starS starS+2*nedge*ngl];
    circS = [circS circS+2*nedge*ngl];
  end
  
  end
  
  function A=SchurBana(P,PW,K,A,starL,circL,starS,circS)
  % use matrix block inversion formula to recursively compute the
  % preconditioner R.
  % 
  % inputs: 
  % P, PW - nontrivial part (i.e., other than the identity matrices) 
  %         prolongation and weighted prolongation matrices
  % K - the system matrix on a type-b mesh along each edge
  %     the type-b mesh contains three chunks, with two small chunks
  %     close to the corner and one big chunk (twice of the small chunk
  %     in the parameter space) away from the corner.
  %
  % A - R_{i-1}, the (i-1)th iteration of R, R is on a type-c mesh
  %     the type-c mesh contains two equal chunks, with the "bad" chunk
  %     close to the corner and the "good" chunk away from the corner.
  %
  %     R is basically a compression from the type-b mesh to the type-c
  %     mesh. When the recursion is done, R on the chunk closest to the
  %     corner contains all the information about the corner singularity in
  %     the sense that 1. during the solve stage, the system matrix is
  %     replace by R on that part of the curve; 2. during the eval stage,
  %     the transformed density rhotilde = R*rho can be integrated with any
  %     other smooth functions using the smooth quadrature rule, i.e., the
  %     Gauss-Legendre rule.
  %
  % starL, circL - indices of "bad" and "good" parts of the system matrix 
  %                bad part - two small chunks close to the corner
  %                good part - one big chunk away from the corner
  %
  % starS, circS - indices of bad and good parts of the preconditioner R
  %                bad part - the chunk close to the corner
  %                good part - the chunk away to the corner
  %
  % output:
  % A - R_{i}, the ith iteration of R.
  %
  VA=K(circL,starL)*A;
  PTA=PW'*A;
  PTAU=PTA*K(starL,circL);
  DVAUI=inv(K(circL,circL)-VA*K(starL,circL));
  DVAUIVAP=DVAUI*(VA*P);
  A(starS,starS)=PTA*P+PTAU*DVAUIVAP;
  A(circS,circS)=DVAUI;
  A(circS,starS)=-DVAUIVAP;
  A(starS,circS)=-PTAU*DVAUI;
  
  end

  function [R]=Rcomp(ngl,nedge,ndim,Pbc,PWbc,nsub,...
    starL,circL,starS,circS,...
    h0,tablocal,clist,isstart,fcurve,rparslocal)
  % the forward recursion for computing the preconditioner R 
  % in the RCIP method
  
  % size of the system matrix
  nsys = 3*ngl*nedge*ndim;
  
  % size of the preconditioner R
  nR = 2*ngl*nedge*ndim;
  
  ts = zeros(4,nedge);
  chnkrlocal(1,nedge) = chunker();
  for level=1:nsub
    h = h0/2^(nsub-level);

    for i=1:nedge
      if isstart(i)
        ts(:,i) = tablocal(1,i)+[0, 0.5*h(i), h(i), 2*h(i)]';
      else
        ts(:,i) = tablocal(2,i)-[2*h(i),h(i),0.5*h(i),0]';
      end
    end
    % construct local chunks around the corner
    for i=1:nedge
      chnkrlocal(i) = chunkerfuncRCIPlocal(fcurve{clist(i)},ts(:,i),pref);
    end
    % construct the system matrix for local chunks
    [MAT,~,~,~] = chunkerbuildclmmat(chnkrlocal,rparslocal);
    %
    MAT = eye(nsys) + MAT;
    if level==1    %  Dumb and lazy initializer for R, for now
      R=eye(nR);  
    end
    R=SchurBana(Pbc,PWbc,MAT,R,starL,circL,starS,circS);   
  end
  
  end


  starL=17:80;
  circL=[1:16 81:96];
  R0=initializeR(theta,Pbc,PWbc,T,W,starL,circL,starS,circS);
  
  omega=logspace(0,3,1000);
  nplot=length(omega);
  myerr=zeros(nplot,1);
  iter=zeros(nplot,1);
  for k=1:nplot
    npan=round(0.6*omega(k)+18);
    sinter=linspace(0,1,npan+1)';
    sinterdiff=ones(npan,1)/npan;
    [z,zp,zpp,w,wzp,awzp,np]=zinit(theta,sinter,sinterdiff,T,W,npan);
%
% *** The K_coa^\circ matrix is set up ***
    LogC=LogCinit(sinterdiff,T,W,npan,np,1);
    K=Koperinit(LogC,z,zp,zpp,wzp,w,omega(k),np);
    S=Soperinit(LogC,z,zp,awzp,sinterdiff,omega(k),np);
    clear LogC
    Kcirc=K-0.5i*omega(k)*S;
    clear K S
    starind=[np-31:np 1:32];
    Kcirc(starind,starind)=zeros(64);
%
% *** Recursion for the R matrix ***
    R=speye(np);
    R(starind,starind)=Rcomp(theta,omega(k),T,W,Pbc,PWbc,starL,circL,R0, ...
			     LogCloc,nsub,npan);
%
% *** Solving main linear system ***
    rhs=2*besselh(0,omega(k)*abs(zsource-z));
    [rhotilde,iter(k)]=myGMRESR(Kcirc,R,rhs,np,itmax,eps);
    clear Kcirc
%
% *** Post processing ***
    rhohat=R*rhotilde;
    tmp=omega(k)*abs(z-ztarg);
    ufield=rhohat.'*(Ktarg(tmp,ztarg,z,wzp)-0.5i*omega(k)*Starg(tmp,awzp))
    ufiref=besselh(0,omega(k)*abs(zsource-ztarg))
    myerr(k)=abs(ufiref-ufield);
    myplot(omega,myerr,iter,k)
  end    

  function R=Rcomp(theta,omega,T,W,Pbc,PWbc,starL,circL,R,LogC,nsub,npan)
  for level=1:nsub
    [z,zp,zpp,w,wzp,awzp,sidi]=zlocinit(theta,T,W,nsub,level,npan);
    K=Koperinit(LogC,z,zp,zpp,wzp,w,omega,96);
    S=Soperinit(LogC,z,zp,awzp,sidi,omega,96);
    MAT=eye(96)+K-0.5i*omega*S;
    R=SchurBana(Pbc,PWbc,MAT,R,starL,circL);
  end
  
  function M1=MRinit(z,zp,zpp,w,wzp,N)
% *** double layer potential ***   
  M1=zeros(N);
  for m=1:N
    M1(:,m)=imag(wzp(m)./(z(m)-z));
  end
  M1(1:N+1:N^2)=w.*imag(zpp./zp)/2;
  M1=M1/pi;

  function S=Soperinit(LogC,z,zp,awzp,sinterdiff,omega,N)
  S=zeros(N);
  for m=1:N
    S(:,m)=1i/2*besselh(0,omega*abs(z-z(m)))*awzp(m);
  end
  myind=find(LogC);
  S(myind)=S(myind)-2*imag(S(myind)).*LogC(myind)/pi;
  tmp=1i/2-(log(omega/2)+0.5772156649015328606)/pi;
  for m=1:N
    dsd=sinterdiff(fix((m-1)/16)+1)/2;
    S(m,m)=(tmp-(LogC(m,m)+log(dsd*abs(zp(m))))/pi)*awzp(m);
  end

  function K=Koperinit(LogC,z,zp,zpp,wzp,w,omega,N)
  K=zeros(N);
  for m=1:N
    tmp=omega*abs(z-z(m));
    K(:,m)=1i/2*tmp.*besselh(1,tmp).*imag(wzp(m)./(z-z(m)));
  end
  K(1:N+1:N^2)=-imag(zpp./zp).*w/2/pi;  
  myind=find(LogC);
  K(myind)=K(myind)-2*imag(K(myind)).*LogC(myind)/pi;
  
  function u=Starg(tmp,awzp)
  u=1i/4*besselh(0,tmp).*awzp;
  
  function u=Ktarg(tmp,ztarg,z,wzp)
  u=1i/4*tmp.*besselh(1,tmp).*imag(wzp./(ztarg-z));

  function R=initializeR(theta,Pbc,PWbc,T,W,starL,circL)
  [z,zp,zpp,w,wzp]=zlocinit0(theta,T,W);
  K=-MRinit(z,zp,zpp,w,wzp,96);
  MAT=eye(96)+K;
  R=inv(MAT(starL,starL));
  iter=0;
  myerr=1;
  while myerr>eps
    Rold=R;
    R=SchurBana(Pbc,PWbc,MAT,R,starL,circL);
    myerr=norm(R-Rold,'fro')/norm(R,'fro');
    iter=iter+1;
  end
  disp(['Fixed point iterations = ',num2str(iter)])

  function A=SchurBana(P,PW,K,A,starL,circL)
  starS=17:48;
  circS=[1:16 49:64];
  VA=K(circL,starL)*A;
  PTA=PW'*A;
  PTAU=PTA*K(starL,circL);
  DVAUI=inv(K(circL,circL)-VA*K(starL,circL));
  DVAUIVAP=DVAUI*(VA*P);
  A(starS,starS)=PTA*P+PTAU*DVAUIVAP;
  A(circS,circS)=DVAUI;
  A(circS,starS)=-DVAUIVAP;
  A(starS,circS)=-PTAU*DVAUI;
  
  function [x,it]=myGMRESR(A,R,b,n,m,tol)
% *** GMRES with low-threshold stagnation control ***
  V=zeros(n,m+1);
  H=zeros(m);
  cs=zeros(m,1);
  sn=zeros(m,1);
  bnrm2=norm(b);
  V(:,1)=b/bnrm2;
  s=bnrm2*eye(m+1,1);
  for it=1:m                                  
    it1=it+1;                                   
    w=A*(R*V(:,it));
    for k=1:it
      H(k,it)=V(:,k)'*w;
      w=w-H(k,it)*V(:,k);
    end
    H(it,it)=H(it,it)+1;
    wnrm2=norm(w);
    V(:,it1)=w/wnrm2;
    for k=1:it-1                                
      temp     = cs(k)*H(k,it)+conj(sn(k))*H(k+1,it);
      H(k+1,it)=-sn(k)*H(k,it)+cs(k)*H(k+1,it);
      H(k,it)  = temp;
    end
    [cs(it),sn(it)]=rotmatc(H(it,it),wnrm2);     
    H(it,it)= cs(it)*H(it,it)+conj(sn(it))*wnrm2;
    s(it1) =-sn(it)*s(it);                      
    s(it)  = cs(it)*s(it);                         
    myerr=abs(s(it1))/bnrm2;
    if (myerr<=tol)||(it==m)                     
      disp(['predicted residual = ' num2str(myerr)])
      y=triu(H(1:it,1:it))\s(1:it);             
      x=fliplr(V(:,1:it))*flipud(y);
      trueres=norm(x+A*(R*x)-b)/bnrm2;
      disp(['true residual      = ',num2str(trueres)])
      break
    end
  end

  function [c,s]=rotmatc(a,b)
  if  a==0
    c=0;
    s=1;
  else
    temp=b/a;
    c=1/sqrt(1+abs(temp)^2);
    s=temp*c;
  end
  
  function [z,zp,zpp,w,wzp,awzp,np]=zinit(theta,sinter,sinterdiff,T,W,npan)
  np=16*npan;
  s=zeros(np,1);
  w=zeros(np,1);
  for k=1:npan
    myind=(k-1)*16+1:k*16;
    sdif=sinterdiff(k)/2;
    s(myind)=(sinter(k)+sinter(k+1))/2+sdif*T;
    w(myind)=W*sdif;
  end
  z=zfunc(s,theta) ;
  zp=zpfunc(s,theta);
  zpp=zppfunc(s,theta);
% *** some extra presicion gained from symmetry ***
  z(np/2+1:np)=conj(flipud(z(1:np/2)));
  zp(np/2+1:np)=-conj(flipud(zp(1:np/2)));
  zpp(np/2+1:np)=conj(flipud(zpp(1:np/2)));
% *************************************************
  wzp=w.*zp;
  awzp=abs(w.*zp);

  function [z,zp,zpp,w,wzp,awzp,sidi]=zlocinit(theta,T,W,nsub,level,npan)
  denom=2^(nsub-level)*npan;
  s=[T/4+0.25;T/4+0.75;T/2+1.5]/denom;
  w=[W/4;W/4;W/2]/denom;
  w=[flipud(w);w];
  sidi=[1;0.5;0.5;0.5;0.5;1]/denom;
  z=zfunc(s,theta);
  z=[conj(flipud(z));z];
  zp=zpfunc(s,theta);
  zp=[-conj(flipud(zp));zp];
  zpp=zppfunc(s,theta);
  zpp=[conj(flipud(zpp));zpp];
  wzp=w.*zp;
  awzp=abs(w.*zp);

  function [z,zp,zpp,w,wzp]=zlocinit0(theta,T,W)
  s=[T/4+0.25;T/4+0.75;T/2+1.5];
  z=s*exp(-1i*theta/2);	    
  z=[conj(flipud(z));z];	    
  zp=ones(48,1)*exp(-1i*theta/2);
  zp=[-conj(flipud(zp));zp];
  zpp=zeros(48,1);
  zpp=[conj(flipud(zpp));zpp];
  w=[W/4;W/4;W/2];
  w=[flipud(w);w];
  wzp=w.*zp;
  


  function [IP,IPW]=IPinit(T,W)
  % construct the prolongation matrix IP that maps function values
  % on n_{gl} Gauss-Legendre nodes on [-1,1] to function values at 
  % 2n_{gl} Gauss-Legendre, with shifted and scaled n_{gl}
  % Gauss-Legendre nodes on each subinterval [-1,0], [0,1], respectively.
  %
  % IPW is the weighted prolongation matrix acted on the left side. 
  ngl = length(T);
  A=ones(ngl);
  AA=ones(2*ngl,ngl);
  T2=[T-1;T+1]/2;
  W2=[W;W]/2;
  for k=2:ngl
    A(:,k)=A(:,k-1).*T;
    AA(:,k)=AA(:,k-1).*T2;   
  end
  IP=AA/A;
  IPW=IP.*(W2*(1./W)');
  
  function Pbc=Pbcinit(IP,nedge,ndim)
  % construct the nontrivial part of the prolongation matrix for the whole
  % system. Suppose that nedge is the number of edges meeting at the 
  % corner, ndim is the number of equations, n_{gl} is the number of 
  % Gauss-Legendre nodes on each chunk. Then Pbc is a block diagonal 
  % matrix with nedge*ndim diagonal blocks, where each block is of the size
  % 2n_{gl} \times n_{gl}
  Pbc=kron(eye(nedge*ndim),IP);
  
  function T=Tinit16
% *** 16-point Gauss-Legendre nodes ***  
  T=zeros(16,1);
  T( 1)=-0.989400934991649932596154173450332627;
  T( 2)=-0.944575023073232576077988415534608345;
  T( 3)=-0.865631202387831743880467897712393132;
  T( 4)=-0.755404408355003033895101194847442268;
  T( 5)=-0.617876244402643748446671764048791019;
  T( 6)=-0.458016777657227386342419442983577574;
  T( 7)=-0.281603550779258913230460501460496106;
  T( 8)=-0.095012509837637440185319335424958063;
  T( 9)= 0.095012509837637440185319335424958063;
  T(10)= 0.281603550779258913230460501460496106;
  T(11)= 0.458016777657227386342419442983577574;
  T(12)= 0.617876244402643748446671764048791019;
  T(13)= 0.755404408355003033895101194847442268;
  T(14)= 0.865631202387831743880467897712393132;
  T(15)= 0.944575023073232576077988415534608345;
  T(16)= 0.989400934991649932596154173450332627;

  function W=Winit16
% *** 16-point Gauss-Legendre weights ***  
  W=zeros(16,1); 
  W( 1)= 0.027152459411754094851780572456018104;
  W( 2)= 0.062253523938647892862843836994377694;
  W( 3)= 0.095158511682492784809925107602246226;
  W( 4)= 0.124628971255533872052476282192016420;
  W( 5)= 0.149595988816576732081501730547478549;
  W( 6)= 0.169156519395002538189312079030359962;
  W( 7)= 0.182603415044923588866763667969219939;
  W( 8)= 0.189450610455068496285396723208283105;
  W( 9)= 0.189450610455068496285396723208283105;
  W(10)= 0.182603415044923588866763667969219939;
  W(11)= 0.169156519395002538189312079030359962;
  W(12)= 0.149595988816576732081501730547478549;
  W(13)= 0.124628971255533872052476282192016420;
  W(14)= 0.095158511682492784809925107602246226;
  W(15)= 0.062253523938647892862843836994377694;
  W(16)= 0.027152459411754094851780572456018104;
  