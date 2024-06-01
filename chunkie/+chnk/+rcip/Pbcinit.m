function Pbc=Pbcinit(IP,nedge,ndim)
  %CHNK.RCIP.Pbcinit
  %
  % construct the nontrivial part of the prolongation matrix for the whole
  % system. Suppose that nedge is the number of edges meeting at the 
  % corner, ndim is the number of equations, n_{gl} is the number of 
  % Gauss-Legendre nodes on each chunk. Then Pbc is a block diagonal 
  % matrix with nedge diagonal blocks, where each block is of the size
  % 2n_{gl} \times n_{gl}
  
  % Pbc=kron(eye(nedge*ndim),IP); % the other order - one big block after
  % another

  % author: Johan Helsing (part of the RCIP tutorial)
  Pbc=kron(eye(nedge),kron(IP,eye(ndim)));
