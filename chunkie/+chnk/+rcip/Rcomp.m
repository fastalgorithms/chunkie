function [R]=Rcomp(ngl,nedge,ndim,Pbc,PWbc,nsub,...
    starL,circL,starS,circS,ilist,...
    h0,isstart,fcurve,glxs,fkern)
  % carry out the forward recursion for computing the preconditioner R 
  % in the RCIP method
  %
  % A more general version of Shidong's rcip version
  %
  % Function is passed as a handle, number of equations is given by
  % ndim
  %
  % Kernel on input takes in arguments (chnkrlocal,ilistl);
  % 
  % Note that matrix must be scaled to have identity on the diagonal,
  % will not work with scaled version of identity
  
  
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
        ts(:,i) =  [0, 0.5, 1, 2]*h(i);
      else
        ts(:,i) = -[2, 1, 0.5, 0]*h(i);
      end
    end
    % construct local chunks around the corner
    for i=1:nedge
      chnkrlocal(i) = chnk.rcip.chunkerfunclocal(fcurve{i},ts(:,i),pref,glxs);
    end
%     figure
%     clf
%     plot(chnkrlocal,'r.')
%     axis equal
%     pause
    % construct the system matrix for local chunks
    if level == 1
      ilistl = [];
    else
      ilistl = ilist;
    end
    
    opts = [];
    % test for opdims ~= [1,1]
    MAT = chunkermat(chnkrlocal,fkern,opts,ilistl);
    
    %
    MAT = eye(nsys) + MAT;
    if level==1    %  Dumb and lazy initializer for R, for now
      %R=eye(nR); 
      R = inv(MAT(starL,starL));
    end
    R=chnk.rcip.SchurBana(Pbc,PWbc,MAT,R,starL,circL,starS,circS);   
  end
  
  end
