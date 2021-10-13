function [clmparams] = setup()
clmparams = [];

ndomain = 5; % number of domains
ncurve = 10; % number of curve segments

rn = zeros(ndomain,1); 
% rn(i) is the index of refraction of the ith domain
rn(1) = 1.0;
rn(2) = 1.4;
rn(3) = 1.5;
rn(4) = 1.6;
rn(5) = 1.55;

% lambda0 is the wavelength of the incident light in vacuum
lambda0 = 0.55/4; % green light wavelength in nm

% k(i) is the wave number for the ith domain
k = rn/lambda0;

% coefficients in the boundary conditions on normal derivatives
coef = 1./rn.^2;

% domain indices for each curve
c = zeros(2,ncurve);
% interior domain for the ith curve
c(1,1:2) = 1; c(1,3:4) = 3; c(1,[5, 6, 7]) = 4; c(1, [8, 9, 10]) = 5; 
% exterior domain for the ith curve 
c(2,1:2) = 2; c(2,3) = 1; c(2,4) = 4; c(2,[5,6,8,9,10]) = 2; c(2,7) = 5;


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

% two circular arcs for the center eye for now
theta = zeros(1,3);
% upper curve opening angle
theta(1) = pi/2.4;

% parameters for the complexification of left and right flat parts.
c1 = log(1/eps)/min(k(1:2));
c2 = c1/2;

% curve parameters
n0 = 50; % average number of wavelength of each curve 
a = -lambda0*n0;
b = -a;

% length for complexification
C = 8*c2;

% curve #1 -- left flat curve [-L(1),a] 
L(1) = C-a;

% curve #2 -- right flat curve [b,L(2)]
L(2) = b+C;

% curve #3 -- circular arc on the real x-axis [a,b]

% curve #4 -- sine curve on [a,b]
n4 = 3; A4 = 0.5;

% curve #5 - #9 -- line segments, specified by vertices
vert = zeros(2,6);
vert(:,1) = [a;0];
vert(:,2) = [b;0];

L(3) = lambda0*n0/2;
theta(2) = pi/8;

% curve #5 -- connect vert(:,1) and vert(:,3)
vert(1,3) = vert(1,1) + L(3)*cos(3*pi/2+theta(2));
vert(2,3) = vert(2,1) + L(3)*sin(3*pi/2+theta(2));

% curve #6 -- connect vert(:,4) and vert(:,2)
vert(1,4) = vert(1,2) + L(3)*cos(3*pi/2-theta(2));
vert(2,4) = vert(2,2) + L(3)*sin(3*pi/2-theta(2));

% curve #7 -- connect vert(:,3) and vert(:,4)
L(4) = lambda0*n0;
% curve #8 -- connect vert(:,3) and vert(:,5)
vert(1,5) = vert(1,3) + L(4)*cos(3*pi/2+theta(2));
vert(2,5) = vert(2,3) + L(4)*sin(3*pi/2+theta(2));

% curve #9 -- connect vert(:,6) and vert(:,4)
vert(1,6) = vert(1,4) + L(4)*cos(3*pi/2-theta(2));
vert(2,6) = vert(2,4) + L(4)*sin(3*pi/2-theta(2));

% curve #10 -- circular arc connecting vert(:,5) and vert(:,6)
theta(3) = 3*pi/4;

% discretize the boundary
tab = zeros(2,ncurve);
tab(:,1) = [-L(1);a];
tab(:,2) = [b;L(2)];

tab(:,3) = [-theta(1)/2;theta(1)/2];
tab(:,10) = [-theta(3)/2;theta(3)/2];

tab(:,4) = [0;n4*pi];

for i=5:9
  tab(1,i) = 0;
  tab(2,i) = 1;
end

cparams = cell(1,ncurve);

for i=1:ncurve
  cparams{i}.ta = tab(1,i);
  cparams{i}.tb = tab(2,i);
  cparams{i}.ifclosed = false;
end

cpars = cell(1,ncurve);
cpars{1}.L = L(1); cpars{1}.c1=c1; cpars{1}.c2=c2;
cpars{2}.L = L(2); cpars{2}.c1=c1; cpars{2}.c2=c2;

cpars{3}.v0 = vert(:,1); cpars{3}.v1 = vert(:,2); 
cpars{3}.theta = theta(1); cpars{3}.ifconvex = 0;

cpars{4}.a = a; cpars{4}.b = b; cpars{4}.A = A4; cpars{4}.n = n4;

cpars{10}.v0 = vert(:,5); cpars{10}.v1 = vert(:,6); 
cpars{10}.theta = theta(3); cpars{10}.ifconvex = 1;

cpars{5}.v0 = vert(:,1); cpars{5}.v1 = vert(:,3);
cpars{6}.v0 = vert(:,4); cpars{6}.v1 = vert(:,2);
cpars{7}.v0 = vert(:,3); cpars{7}.v1 = vert(:,4);
cpars{8}.v0 = vert(:,3); cpars{8}.v1 = vert(:,5);
cpars{9}.v0 = vert(:,6); cpars{9}.v1 = vert(:,4);

ncorner = 6;
issymmetric = 1;
corners = cell(1,ncorner);

for icorner = 1:ncorner
  corners{icorner} = [];
  if icorner==1
    clist = [1, 3, 4, 5]; % curve ID starting or ending at the corner
    isstart = [0, 0, 1, 1]; % 1 -- corner is the starting point; 0 -- corner is the end point
    nedge = 4; % number of edges
  elseif icorner == 2
    clist = [2, 3, 4, 6];
    isstart = [1, 1, 0, 0];
    nedge = 4; % quadruple junction
  elseif icorner == 3
    clist = [5, 7, 8];
    isstart = [0, 1, 1];
    nedge = 3; % triple junction
  elseif icorner == 4
    clist = [6, 7, 9];
    isstart = [1, 0, 0];
    nedge = 3;
  elseif icorner == 5
    clist = [8, 10];
    isstart = [0, 1];
    nedge = 2;
  elseif icorner == 6
    clist = [9, 10];
    isstart = [1, 0];
    nedge = 2;
  end
  
  corners{icorner}.clist = clist;
  corners{icorner}.isstart = isstart;
  corners{icorner}.nedge = nedge;
end
  
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
