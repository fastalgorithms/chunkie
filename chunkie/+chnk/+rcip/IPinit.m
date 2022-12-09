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