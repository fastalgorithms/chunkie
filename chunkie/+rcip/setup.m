  
  function [Pbc,PWbc,starL,circL,starS,circS] = setup(ngl,ndim,nedge,isstart) 
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
  [IP,IPW]=rcip.IPinit(T,W);
  Pbc = rcip.Pbcinit(IP,nedge,ndim);
  PWbc = rcip.Pbcinit(IPW,nedge,ndim);
  
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
  
  for i=1:ndim-1
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
      starS = [starS indb0+2*(i-1)*ngl];
      circS = [circS indg0+2*(i-1)*ngl]; 
    end
  end
  
  for i=1:ndim-1
    starS = [starS starS+2*nedge*ngl];
    circS = [circS circS+2*nedge*ngl];
  end
  
  end
  
