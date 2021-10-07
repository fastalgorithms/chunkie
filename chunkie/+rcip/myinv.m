  function Mi=myinv(M)
% *** equation (9) in Henderson & Searle ***
  np=length(M)/2;
  A=M(1:np,1:np);
  U=M(1:np,np+1:2*np);
  V=M(np+1:2*np,1:np);
  D=M(np+1:2*np,np+1:2*np);  
  Mi=[inv(A-U/D*V) -A\U/(D-V/A*U);
      -D\V/(A-U/D*V) inv(D-V/A*U)];