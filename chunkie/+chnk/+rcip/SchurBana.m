  function A=SchurBana(P,PW,K,A,starL,circL,starS,circS)
  %CHNK.RCIP.SCHURBANA block inverse used in RCIP
  % 
  % uses a matrix block inversion formula to recursively compute the
  % preconditioner R.
  %
  % R_{i} =   [PW^T 0] [A^(-1)   U]^(-1) [P  0]
  %           [ 0   I] [V        D]      [0  I]  
  %
  %           2N x 3N      3N x 3N      3N x 2N
  %
  % where A = R_{i-1}, U = K(bad,good), V = K(good,bad), D = K(good,good) 
  % and bad refers to the two close panels on a type b mesh and good the 
  % larger far panel. 
  % 
  % This routine uses the Schur-Banachiewicz speedup as suggested by
  % Helsing (see e.g. the RCIP Tutorial eq 31):
  %
  % R_{i} =   [M1 M2]
  %           [M3 M4]
  %
  % M1 = PW^T A P + PW^T A U (D - VAU)^(-1) VAP
  % M2 = -PW^TAU(D-VAU)^(-1) 
  % M3 = -(D-VAU)^(-1)VAP
  % M4 = (D-VAU)^(-1)
  %
  % this reduces a 3Nx3N matrix inverse to a NxN matrix inverse and
  % several mat-mats
  %
  % inputs: 
  % P, PW - nontrivial part (i.e., other than the identity matrices) 
  %         prolongation and weighted prolongation matrices from type c to
  %         type b meshes
  % K - the system matrix on a type-b mesh along each edge.
  %     the type-b mesh contains three chunks, with two small chunks
  %     close to the corner and one big chunk (twice of the small chunk
  %     in the parameter space) away from the corner.
  %
  % A - R_{i-1}, the (i-1)th iteration of R, R is on a type-c mesh.
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
  
  % author: Johan Helsing (part of the RCIP tutorial)
    
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
