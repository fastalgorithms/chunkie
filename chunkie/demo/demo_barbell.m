%DEMO_BARBELL
%
% Define a polygonal "barbell" shaped domain with 
% prescribed constant boundary data on each edge. 
% Solves the corresponding Laplace Dirichlet problem
% using corner rounding and smoothing the boundary data 
% across the rounded corners
%

iseed = 8675309;
rng(iseed);
addpaths_loc();

% define geometry and boundary conditions
% (vertices defined by another function)
verts = chnk.demo.barbell(2.0,2.0,1.0,1.0);
nv = size(verts,2);
edgevals = rand(1,nv);

% parameters for curve rounding/chunking routineime to oversample boundary
cparams = [];
cparams.widths = 0.01*ones(size(verts,2),1);% width to cut about each corner
cparams.eps = 1e-12; % tolerance at which to resolve curve

% parameters for chunk discretization
p.k = 16; p.dim = 2;

% call smoothed-polygon chunking routine
% a smoothed version of edgevals is returned in 
% chnkr.data
chnkr = chunkpoly(verts,cparams,p,edgevals);
chnkr = chnkr.refine(); chnkr = chnkr.sort();
assert(checkadjinfo(chnkr) == 0);

% plot smoothed geometry and data

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal

figure(2)
clf
nchplot = 1:chnkr.nch;
x = chnkr.r(1,:,nchplot);
y =chnkr.r(2,:,nchplot);
z = chnkr.data(1,:,nchplot);
plot3(x(:),y(:),z(:))

% solve and visualize the solution

% build laplace dirichlet matrix

fkern = @(s,t,stau,ttau) chnk.lap2d.kern(s,t,stau,ttau,'D');
opdims(1) = 1; opdims(2) = 1;
intparams.intorder = chnkr.k;
start = tic; D = chunkskernmat(chnkr,fkern,opdims,intparams);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(chnkr.k*chnkr.nch) + D;

rhs = chnkr.data(1,:); rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
nplot = 100;
hx = (rmax(1)-rmin(1))/nplot;
hy = (rmax(2)-rmin(2))/nplot;
xtarg = linspace(rmin(1)+hx/2,rmax(1)-hx/2,nplot); 
ytarg = linspace(rmin(2)+hy/2,rmax(2)-hy/2,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic;
chnkr2 = chnkr;
chnkr2 = chnkr2.makedatarows(1);
chnkr2.data(2,:) = sol(:);
optref = []; optref.nover = 2;
chnkr2 = chnkr2.refine(optref);
sol2 = chnkr2.data(2,:);
t1 = toc(start);

fprintf('%5.2e s : time to oversample boundary\n',t1)

%

start = tic; in = chunkerinflam(chnkr,targets); t1 = toc(start);

fprintf('%5.2e s : time to find points in domain\n',t1)

% compute layer potential based on oversample boundary

wts2 = whts(chnkr2);

matfun = @(i,j) kernbyindexr(i,j,targets(:,in),chnkr2,wts2,fkern,opdims);
[pr,ptau,pw,pin] = proxy_square_pts();

pxyfun = @(rc,rx,cx,slf,nbr,l,ctr) proxyfunr(rc,rx,slf,nbr,l,ctr,chnkr2,wts2, ...
    fkern,opdims,pr,ptau,pw,pin);

xflam = chnkr2.r(:,:);

start = tic; F = ifmm(matfun,targets(:,in),xflam,200,1e-14,pxyfun); 
t1 = toc(start);
fprintf('%5.2e s : time for ifmm form (for plotting)\n',t1)
start = tic;
Dsol = ifmm_mv(F,sol2(:),matfun); t1 = toc(start);
fprintf('%5.2e s : time for ifmm apply (for plotting)\n',t1)

% 

figure(3)
clf
zztarg = nan(size(xxtarg));
zztarg(in) = Dsol;
h=surf(xxtarg,yytarg,zztarg);
set(h,'EdgeColor','none')