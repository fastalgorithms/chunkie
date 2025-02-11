function [Pbc,PWbc,starL,circL,starS,circS,ilist,starL1,circL1] = ...
      setup(ngl,ndim,nedge,isstart) 
  %CHNK.RCIP.SETUP
  %
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
  % starL1, circL1 - bad and good indices for arrays of nodes, normals, etc
  % starS, circS - bad and good indices for the preconditioner R
  %

  % author: Shidong Jiang, parts drawn from Johan Helsing's RCIP tutorial
  % modified: Manas Rachh, Jeremy Hoskins
  
  [T,W] = lege.exps(ngl);  
  [IP,IPW]=chnk.rcip.IPinit(T,W);
  Pbc = chnk.rcip.Pbcinit(IP,nedge,ndim);
  PWbc = chnk.rcip.Pbcinit(IPW,nedge,ndim);
  
  ilist = zeros(2,nedge);
  % starL - bad indices for the system matrix M
  % circL - good indices for the system matrix M
  starL = [];
  circL = [];
  starL1 = [];
  circL1 = [];
  indg1 = 2*ngl*ndim + (1:ngl*ndim);
  indb1 = 1:2*ngl*ndim;
  indg11 = 2*ngl+ (1:ngl);
  indb11 = 1:2*ngl;
  
  indg0 = 1:ngl*ndim;
  indb0 = ngl*ndim + (1:2*ngl*ndim);
  indg01 = 1:ngl;
  indb01 = ngl + (1:2*ngl);
  
  for i=1:nedge
    if isstart(i) 
      starL = [starL indb1+3*(i-1)*ngl*ndim];
      circL = [circL indg1+3*(i-1)*ngl*ndim];
      starL1 = [starL1 indb11+3*(i-1)*ngl];
      circL1 = [circL1 indg11+3*(i-1)*ngl];
      ilist(:,i) = [1, 2];
    else
      starL = [starL indb0+3*(i-1)*ngl*ndim];
      circL = [circL indg0+3*(i-1)*ngl*ndim];
      starL1 = [starL1 indb01+3*(i-1)*ngl];
      circL1 = [circL1 indg01+3*(i-1)*ngl];
      ilist(:,i) = [2, 3];
    end
  end
  
%   for i=1:ndim-1
%     starL = [starL starL+3*nedge*ngl]; % type-b mesh has 3 chunks
%     circL = [circL circL+3*nedge*ngl];
%   end
  
  % starS - bad indices for the preconditioner R
  % circS - good indices for the preconditioner R
  starS = [];
  circS = [];
  indg1 = ngl*ndim + (1:ngl*ndim);
  indb1 = 1:ngl*ndim;
  
  indg0 = 1:ngl*ndim;
  indb0 = ngl*ndim + (1:ngl*ndim);
  
  for i=1:nedge
    if isstart(i)
      starS = [starS indb1+2*(i-1)*ngl*ndim]; % type-c mesh has 2 chunks
      circS = [circS indg1+2*(i-1)*ngl*ndim];
    else
      starS = [starS indb0+2*(i-1)*ngl*ndim];
      circS = [circS indg0+2*(i-1)*ngl*ndim]; 
    end
  end
  
%   for i=1:ndim-1
%     starS = [starS starS+2*nedge*ngl];
%     circS = [circS circS+2*nedge*ngl];
%   end
  
  end
  
