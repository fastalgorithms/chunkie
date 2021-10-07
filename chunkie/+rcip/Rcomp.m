  function [R]=Rcomp(ngl,nedge,ndim,Pbc,PWbc,nsub,...
    starL,circL,starS,circS,...
    h0,tablocal,clist,isstart,fcurve,rparslocal)
  % carry out the forward recursion for computing the preconditioner R 
  % in the RCIP method
  
  pref = []; 
  pref.k = ngl;
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
      chnkrlocal(i) = rcip.chunkerfunclocal(fcurve{clist(i)},ts(:,i),pref);
    end
%     figure
%     clf
%     plot(chnkrlocal,'r.')
%     axis equal
%     pause
    % construct the system matrix for local chunks
    [MAT,~,~,~] = clm.buildmat(chnkrlocal,rparslocal);
    %
    MAT = eye(nsys) + MAT;
    if level==1    %  Dumb and lazy initializer for R, for now
      %R=eye(nR); 
      R = inv(MAT(starL,starL));
    end
    R=rcip.SchurBana(Pbc,PWbc,MAT,R,starL,circL,starS,circS);   
  end
  
  end