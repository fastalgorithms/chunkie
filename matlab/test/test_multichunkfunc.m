
%TESTMULTICHUNKFUNC
%
% This file tests whether or not multichunkfunc is properly wrapped

ncomp = 2;
k = 16;
lpars = 1000; % make long enough for domain description (at least 2*nmodes+2 +20)

type = 'rmodes';

eps = 1.0e-13*ones(ncomp,1);
ifclosed = ones(ncomp,1);
novers = 2*ones(ncomp,1);
chsmall = ones(ncomp,1);
tas = 0*ones(ncomp,1);
tbs = 2*pi*ones(ncomp,1);

pars = zeros(lpars,1);
ipars = zeros(ncomp,1);

ccs = zeros(2,ncomp);

% describe boundary components in ipars and pars

% component 1
nmodes = 6;
cx = 0.0; cy = 0.0;
ccs(1,1) = cx; ccs(2,1) = cy;
ipars(1) = 1; % starting index of component 1
ii = ipars(1);
pars(ii) = nmodes + 0.1;
pars(ii+1) = 6 + 0.1; % always 6 (internal convention)
pars(ii+2) = 1 + 0.1; % sign of orientation
pars(ii+3) = cx;
pars(ii+4) = cy;
pars(ii+5) = 0.25; % coeff of constant zero mode
% cos and sin coeffs for higher modes
pars(ii+6) = 0; pars(ii+7) = 0;
pars(ii+8) = 0; pars(ii+9) = 0;
pars(ii+10) = 0; pars(ii+11) = 0.01;
pars(ii+12) = 0.02; pars(ii+13) = 0;
pars(ii+14) = 0.01; pars(ii+15) = 0;
pars(ii+16) = 0; pars(ii+17) = 0; % nmodes = 6

%component 2
nmodes = 5;
cx = 0.0; cy = 0.0;
ccs(1,2) = cx; ccs(2,2) = cy;
ipars(2) = ii+19;
ii = ipars(2);
pars(ii) = nmodes + 0.1;
pars(ii+1) = 6 + 0.1; % always 6 (internal convention)
pars(ii+2) = 1 + 0.1; % sign of orientation
pars(ii+3) = cx;
pars(ii+4) = cy;
pars(ii+5) = 0.05; % coeff of constant zero mode
% cos and sin coeffs for higher modes
pars(ii+6) = 0; pars(ii+8) = 0;
pars(ii+8) = 0.005; pars(ii+9) = 0;
pars(ii+10) = 0; pars(ii+11) = 0.005;
pars(ii+12) = 0; pars(ii+13) = 0;
pars(ii+14) = 0.005; pars(ii+15) = 0; % nmodes = 5

geopars.eps = eps;
geopars.ifclosed = ifclosed;
geopars.chsmall = chsmall;
geopars.tas = tas;
geopars.tbs = tbs;
geopars.pars = pars;
geopars.ipars = ipars;
geopars.k = k;
geopars.novers = novers;
geopars.type = type;
geopars.ncomp = ncomp;

[wgeos,info] = multichunkfunc(geopars);
%wgeo = multichunkgetcomp(wgeos,1);

[wgeo,nchs] = multichunkmergepack(wgeos);

chunker = chunkunpack(wgeo);

xs = chunker.chunks(1,:,:); xs = xs(:);
ys = chunker.chunks(2,:,:); ys = ys(:);

figure(1)
clf
scatter(xs,ys)