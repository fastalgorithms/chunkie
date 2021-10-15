function [clmparams] = setup2()
clmparams = [];

ndomain = 2; % number of domains
ncurve = 1; % number of curve segments

rn = zeros(ndomain,1); 
% rn(i) is the index of refraction of the ith domain
rn(1) = 4.5;
rn(2) = 4.2+1i*1e-1;

% k0 is the wave number in vacuum
k0 = 1;

% k(i) is the wave number for the ith domain
k = k0*rn;


% coefficients in the boundary conditions on normal derivatives
%coef = rn.^2; % TM polarization
coef = ones(size(rn)); % TE polarization

% domain indices for each curve
c = zeros(2,ncurve);
% interior domain for the ith curve
c(1) = 1; 
% exterior domain for the ith curve 
c(2) = 2;

k1 = k(1); % wave numbers for the interior domain
k2 = k(2); % wave numbers for the exterior domain

% parameters for the complexification of left and right flat parts.
% curve parameters
L = 10;
c1 = 8;
c2 = 3;

% discretize the boundary
tab = zeros(2,ncurve);
tab(:,1) = [-L;L];

cparams = cell(1,ncurve);

for i=1:ncurve
  cparams{i}.ta = tab(1,i);
  cparams{i}.tb = tab(2,i);
  cparams{i}.ifclosed = false;
end


cpars = cell(1,ncurve);
cpars{1}.L = L; cpars{1}.c1=c1; cpars{1}.c2=c2;

clist = cell(1,ndomain);
for i=1:ndomain
  clist{i} = [];
  for j=1:ncurve
    if c(1,j)==i
      clist{i} = [clist{i} j];
    elseif c(2,j)==i
      clist{i} = [clist{i} -j];
    end
  end
end
clmparams.clist = clist;


% number of chunks on each curve
% should be proportional to the length*wavenumber
nch = zeros(1,ncurve);
nch(1) = 2*L;

% allocate sources for point source test
src = [1.1 0;3 -2];



clmparams.k = k;
clmparams.c = c;

clmparams.k1 = k1;
clmparams.k2 = k2;
clmparams.coef = coef;
clmparams.cpars = cpars;
clmparams.cparams = cparams;
clmparams.ncurve = ncurve;
clmparams.ndomain = ndomain;


clmparams.nch = nch;
clmparams.src = src;
