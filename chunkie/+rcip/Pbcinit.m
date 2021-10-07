  function Pbc=Pbcinit(IP,nedge,ndim)
  % construct the nontrivial part of the prolongation matrix for the whole
  % system. Suppose that nedge is the number of edges meeting at the 
  % corner, ndim is the number of equations, n_{gl} is the number of 
  % Gauss-Legendre nodes on each chunk. Then Pbc is a block diagonal 
  % matrix with nedge*ndim diagonal blocks, where each block is of the size
  % 2n_{gl} \times n_{gl}
  Pbc=kron(eye(nedge*ndim),IP);