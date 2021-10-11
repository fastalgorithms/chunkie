function [M,np,alpha1,alpha2] = buildmat_fast(chnkr,rpars,opts,opdims,...
  glwts,iglist,logquad)
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

% diagonal constant for each curve
alpha1 = zeros(1,ncurve);
alpha2 = zeros(1,ncurve);
for i=1:ncurve
  alpha1(i) = 2/(coef(c(1,i))+coef(c(2,i)));
  alpha2(i) = 2/(1/coef(c(1,i))+1/coef(c(2,i)));
end

quad = [];
if isfield(opts, 'quad')
  quad = opts.quad;
end

for i=1:ncurve % target curve id
  c1 = coef(c(1,i));
  c2 = coef(c(2,i));
  %
  % define kernels
  allk1 =  @(s,t) chnk.helm2d.kern(k1(i),s,t,'all',1);
  allk2 =  @(s,t) chnk.helm2d.kern(k2(i),s,t,'all',1);

  indi1 = sum(nch(1:i-1))*2*ngl+(1:2:2*nch(i)*ngl);
  indi2 = sum(nch(1:i-1))*2*ngl+(2:2:2*nch(i)*ngl);
  
  ni1 = 1:2:2*nch(i)*ngl;
  ni2 = 2:2:2*nch(i)*ngl;
  
  for j=1:ncurve % source curve id
    if j==i
      % build matrices for 8 layer potentials
      jlist = [];
      if ~isempty(iglist)
        jlist = iglist(:,j);
      end
      
      logquad.omega = k1(i);
      M1 = chunkermat_fast(chnkr(i),allk1,opts,glwts,jlist,logquad);
     
      logquad.omega = k2(i);
      M2 = chunkermat_fast(chnkr(i),allk2,opts,glwts,jlist,logquad);
      
      M(indi1,indi1) =  alpha1(i)*(c2*M2(ni1,ni1)-c1*M1(ni1,ni1)); % D
      M(indi1,indi2) =  alpha1(i)*(M2(ni1,ni2)-M1(ni1,ni2)); % S
      
      M(indi2,indi1) = -alpha2(i)*(M2(ni2,ni1)-M1(ni2,ni1)); % D'
      M(indi2,indi2) = -alpha2(i)*(1/c2*M2(ni2,ni2)-1/c1*M1(ni2,ni2)); % S'  
    else
      indj1 = sum(nch(1:j-1))*2*ngl+(1:2:2*nch(j)*ngl);
      indj2 = sum(nch(1:j-1))*2*ngl+(2:2:2*nch(j)*ngl);
      
      nj1 = 1:2:2*nch(j)*ngl;
      nj2 = 2:2:2*nch(j)*ngl;

      ilist = [];
      jlist = [];
%       if ~isempty(iglist)
%         ilist = iglist(:,i);
%         jlist = iglist(:,j);
%       end
      
      M1 = chunkermat_smooth(chnkr(j),chnkr(i),allk1,opdims,glwts,jlist,ilist);
      M2 = chunkermat_smooth(chnkr(j),chnkr(i),allk2,opdims,glwts,jlist,ilist);
      
      M(indi1,indj1) =  alpha1(i)*(c2*M2(ni1,nj1)-c1*M1(ni1,nj1));
      M(indi1,indj2) =  alpha1(i)*(M2(ni1,nj2)-M1(ni1,nj2));
      
      M(indi2,indj1) = -alpha2(i)*(M2(ni2,nj1)-M1(ni2,nj1));
      M(indi2,indj2) = -alpha2(i)*(1/c2*M2(ni2,nj2)-1/c1*M1(ni2,nj2));
    end
  end
end

end