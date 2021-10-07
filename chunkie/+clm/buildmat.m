function [M,np,alpha1,alpha2] = buildmat(chnkr,rpars)
% build the system matrix for the interface problem
% inputs: 
% rpars.k - array of wave numbers for each domain
% rpars.c - (2,ncurve) matrix, where c(1,i) returns the domain index for the
%     interior region w.r.t. the ith curve, and c(2,i) returns the domain
%     index for the exterior region w.r.t. the ith curve.
%
% rpars.coef - array of the reciprocal of coefficients for each domain in the 
%        boundary condition for the normal derivative
% chnkr - array of chnkr objects specifying the whole boundary,
%         For constructing the preconditioner R in the RCIP method, chnkr
%         contains the type-b mesh on each edge
% outputs:
% M - the system matrix
%
ncurve = length(chnkr);
nch = zeros(1,ncurve);

k = rpars.k;
c = rpars.c;
coef = rpars.coef;

for i=1:ncurve
    nch(i) = chnkr(i).nch;
end

k1 = zeros(1,ncurve); % wave numbers for the interior domain
k2 = zeros(1,ncurve); % wave numbers for the exterior domain

for i=1:ncurve
  k1(i) = k(c(1,i));
  k2(i) = k(c(2,i));
end

% total number of discretization points
ngl = chnkr(1).k;
np = sum(nch(1:ncurve))*ngl;

%[~,~,info] = sortinfo(chnkr);
%assert(info.ier == 0);


% now build the system matrix
M = zeros(2*np);

opts = [];

% diagonal constant for each curve
alpha1 = zeros(1,ncurve);
alpha2 = zeros(1,ncurve);
for i=1:ncurve
  alpha1(i) = 2/(coef(c(1,i))+coef(c(2,i)));
  alpha2(i) = 2/(1/coef(c(1,i))+1/coef(c(2,i)));
end

for i=1:ncurve % target curve id
  c1 = coef(c(1,i));
  c2 = coef(c(2,i));
  %
  % define kernels
  sk1 =  @(s,t) chnk.helm2d.kern(k1(i),s,t,'s',1);
  sk2 =  @(s,t) chnk.helm2d.kern(k2(i),s,t,'s',1);

  dk1 =  @(s,t) chnk.helm2d.kern(k1(i),s,t,'d',1);
  dk2 =  @(s,t) chnk.helm2d.kern(k2(i),s,t,'d',1);

  spk1 = @(s,t) chnk.helm2d.kern(k1(i),s,t,'sprime',1);
  spk2 = @(s,t) chnk.helm2d.kern(k2(i),s,t,'sprime',1);

  dpk1 = @(s,t) chnk.helm2d.kern(k1(i),s,t,'dprime',1);
  dpk2 = @(s,t) chnk.helm2d.kern(k2(i),s,t,'dprime',1);
  
  indi1 = sum(nch(1:i-1))*ngl+(1:nch(i)*ngl);
  indi2 = (indi1)+np;
  
  for j=1:ncurve % source curve id
    if j==i
      % build matrices for 8 layer potentials
      S1  = chunkermat(chnkr(i),sk1,opts);
      S2  = chunkermat(chnkr(i),sk2,opts);

      D1  = chunkermat(chnkr(i),dk1,opts);
      D2  = chunkermat(chnkr(i),dk2,opts);

      S1p = chunkermat(chnkr(i),spk1,opts);
      S2p = chunkermat(chnkr(i),spk2,opts);

      D1p = chunkermat(chnkr(i),dpk1,opts);
      D2p = chunkermat(chnkr(i),dpk2,opts);
      
      M(indi1,indi1) =  alpha1(i)*(c2*D2-c1*D1);
      M(indi1,indi2) =  alpha1(i)*(S2-S1);
      
      M(indi2,indi1) = -alpha2(i)*(D2p-D1p);
      M(indi2,indi2) = -alpha2(i)*(1/c2*S2p-1/c1*S1p);
    else
      indj1 = sum(nch(1:j-1))*ngl+(1:nch(j)*ngl);
      indj2 = (indj1)+np;
        
      S1  = chunkermat_smooth(chnkr(j),chnkr(i),sk1);
      S2  = chunkermat_smooth(chnkr(j),chnkr(i),sk2);
      
      D1  = chunkermat_smooth(chnkr(j),chnkr(i),dk1);
      D2  = chunkermat_smooth(chnkr(j),chnkr(i),dk2);
      
      S1p = chunkermat_smooth(chnkr(j),chnkr(i),spk1);
      S2p = chunkermat_smooth(chnkr(j),chnkr(i),spk2);
      
      D1p = chunkermat_smooth(chnkr(j),chnkr(i),dpk1);       
      D2p = chunkermat_smooth(chnkr(j),chnkr(i),dpk2); 
      
      M(indi1,indj1) =  alpha1(i)*(c2*D2-c1*D1);
      M(indi1,indj2) =  alpha1(i)*(S2-S1);
      
      M(indi2,indj1) = -alpha2(i)*(D2p-D1p);
      M(indi2,indj2) = -alpha2(i)*(1/c2*S2p-1/c1*S1p);
      
      
      if 1==2
      % now fix the near interactions between adjacent chunks on different
      % curves
      srcid=[];
      if i==1 && j==4
        targid = nch(i);
        srcid = 1;
      elseif i==4 && j==1
        targid = 1;
        srcid = nch(j);
      elseif (i==1 && j==3) || (i==3 && j==1)
        srcid = nch(j);
        targid = nch(i);
      elseif i==2 && j==4
        srcid = nch(j);
        targid = 1;
      elseif i==4 && j==2
        srcid = 1;
        targid = nch(i);
      elseif (i==2 && j==3) || (i==3 && j==2)
        srcid = 1;
        targid = 1;
      elseif (i==3 && j==4) || (i==4 && j==3)
        srcid = [1, nch(j)];
        targid = [nch(i), 1];
      end
      
      for m=1:length(srcid)
      S1  = chunkermat_targ(chnkr(j),chnkr(i),srcid(m),targid(m),sk1,opdims,'nearlog');
      S2  = chunkermat_targ(chnkr(j),chnkr(i),srcid(m),targid(m),sk2,opdims,'nearlog');
      D1  = chunkermat_targ(chnkr(j),chnkr(i),srcid(m),targid(m),dk1,opdims,'nearlog');
      D2  = chunkermat_targ(chnkr(j),chnkr(i),srcid(m),targid(m),dk2,opdims,'nearlog');      
      S1p = chunkermat_targ(chnkr(j),chnkr(i),srcid(m),targid(m),spk1,opdims,'nearlog');
      S2p = chunkermat_targ(chnkr(j),chnkr(i),srcid(m),targid(m),spk2,opdims,'nearlog');
      D1p = chunkermat_targ(chnkr(j),chnkr(i),srcid(m),targid(m),dpk1,opdims,'nearlog');
      D2p = chunkermat_targ(chnkr(j),chnkr(i),srcid(m),targid(m),dpk2,opdims,'nearlog'); 
      
      idi1 = (sum(nch(1:i-1))+targid(m)-1)*ngl+(1:ngl);
      idi2 = idi1 + np;
      
      idj1 = (sum(nch(1:j-1))+srcid(m)-1)*ngl+(1:ngl);
      idj2 = idj1 + np;
      
      M(idi1,idj1) =  alpha1(i)*(c2*D2-c1*D1);
      M(idi1,idj2) =  alpha1(i)*(S2-S1);
      
      M(idi2,idj1) = -alpha2(i)*(D2p-D1p);
      M(idi2,idj2) = -alpha2(i)*(1/c2*S2p-1/c1*S1p);   
      end
      end
    end
  end
end

end