function [clmparams] = setup4()
clmparams = [];

ndomain = 3; % number of domains
ncurve = 4; % number of curve segments

rn = zeros(ndomain,1); 
% rn(i) is the index of refraction of the ith domain
rn(1) = 1.0;
rn(2) = 1.6+1i*0e-1;
rn(3) = 1.4;

% k0 is the wave number in vacuum
k0 = 4;

% k(i) is the wave number for the ith domain
k = k0*rn;


% coefficients in the boundary conditions on normal derivatives
%coef = rn.^2; % TM polarization
coef = ones(size(rn)); % TE polarization

% domain indices for each curve
c = zeros(2,ncurve);
% interior domain for the ith curve
c(1,1:2) = 1; c(1,3:4) = 3;
% exterior domain for the ith curve 
c(2,1:2) = 2; c(2,3) = 1; c(2,4) = 2;

k1 = zeros(1,ncurve); % wave numbers for the interior domain
k2 = zeros(1,ncurve); % wave numbers for the exterior domain

for i=1:ncurve
  k1(i) = k(c(1,i));
  k2(i) = k(c(2,i));
end

% two circular arcs for the center eye for now
theta = zeros(1,ncurve);
% upper curve opening angle
theta(3) = pi/2.4;
% lower curve opening angle
theta(4) = pi/2.2;

% parameters for the complexification of left and right flat parts.
c1 = log(1d-2/eps)/min(abs(k(1:2)))
c2 = c1/2.5
% curve parameters
% center eye range on the real x-axis [a,b]
a = -1d1;
b = 1d1;
% length for complexification
C = 6*c2
% left flat curve [-L(1),a]
L(1) = C-a;
% right flat curve [b,L(2)]
L(2) = b+C;


% discretize the boundary
tab = zeros(2,ncurve);
tab(:,1) = [-L(1);a];
tab(:,2) = [b;L(2)];
for i=3:4
  tab(:,i) = [-theta(i)/2;theta(i)/2];
end

cparams = cell(1,ncurve);

for i=1:ncurve
  cparams{i}.ta = tab(1,i);
  cparams{i}.tb = tab(2,i);
  cparams{i}.ifclosed = false;
end

vert = zeros(2,2);
vert(:,1) = [a;0];
vert(:,2) = [b;0];

cpars = cell(1,ncurve);
cpars{1}.L = L(1); cpars{1}.c1=c1; cpars{1}.c2=c2;
cpars{2}.L = L(2); cpars{2}.c1=c1; cpars{2}.c2=c2;

cpars{3}.v0 = vert(:,1); cpars{3}.v1 = vert(:,2); 
cpars{3}.theta = theta(3); cpars{3}.ifconvex = 0;
cpars{4}.v0 = vert(:,1); cpars{4}.v1 = vert(:,2); 
cpars{4}.theta = theta(4); cpars{4}.ifconvex = 1;


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

k1 = zeros(1,ncurve); % wave numbers for the interior domain
k2 = zeros(1,ncurve); % wave numbers for the exterior domain

for i=1:ncurve
  k1(i) = k(c(1,i));
  k2(i) = k(c(2,i));
end


ncorner = 2;
issymmetric = 1;
corners = cell(1,ncorner);

for icorner = 1:ncorner
  corners{icorner} = [];
  if icorner==1
    clist = [1, 3, 4];
    isstart = [0, 0, 1];
    nedge = 3;
  else
    clist = [2, 3, 4];
    isstart = [1, 1, 0];
    nedge = 3;
  end
  corners{icorner}.clist = clist;
  corners{icorner}.isstart = isstart;
  corners{icorner}.nedge = nedge;
end
  

% number of chunks on each curve
% should be proportional to the length*wavenumber
nch = zeros(1,ncurve);

n0 = 6;
fac = 1.2;

lambda = 2*pi/max(abs(k(1)),abs(k(2)));
nch(1) = round(fac*(a+L(1))/lambda) + n0;
nch(2) = round(fac*(L(2)-b)/lambda) + n0;

n0 = 16;
lambda = 2*pi/max(abs(k(1)),abs(k(3)));
nch(3) = round(fac*(b-a)/2/sin(theta(3)/2)/lambda) + n0;
lambda = 2*pi/max(abs(k(2)),abs(k(3)));
nch(4) = round(fac*(b-a)/2/sin(theta(4)/2)/lambda) + n0;

% allocate sources for point source test
src = [-1, 1.3, 0; 4.4,-4.8,0];



clmparams.k = k;
clmparams.c = c;

clmparams.k1 = k1;
clmparams.k2 = k2;
clmparams.coef = coef;
clmparams.cpars = cpars;
clmparams.cparams = cparams;
clmparams.ncurve = ncurve;
clmparams.ndomain = ndomain;
clmparams.ncorner = ncorner;
clmparams.issymmetric = issymmetric;
clmparams.corners = corners;
clmparams.vert = vert;

clmparams.nch = nch;
clmparams.src = src;
