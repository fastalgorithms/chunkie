%DEMO_SCATTER_SLIT
%
% Define a scattering problem on a slit-like domain and solve
%

iseed = 8675309;
rng(iseed);
addpaths_loc();

% incident wave and other problem definitions

src0 = [0.7;-5];
zk = 2.3 + 1i*0.001;
strengths = 1.0;

% define geometry

% chunkpoly is based on vertices

% 4 wide, 1 tall rectangle at origin

verts = [ [-2;-0.5],[2;-0.5],[2;0.5],[-2;0.5] ];

cparams = [];
cparams.eps = 1.0e-5;
cparams.rounded = true;
pref = []; 
pref.k = 16;
start = tic; chnkr = chunkerpoly(verts,cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

assert(checkadjinfo(chnkr) == 0);
refopts = []; refopts.maxchunklen = 4.0/abs(zk);
chnkr = chnkr.refine(refopts); chnkr = chnkr.sort();

% make 2 shifted copies of this chnkr and merge them

chnkr1 = chnkr;
chnkr1.r = chnkr1.r - [3;0]; % left rectangle
chnkr2 = chnkr;
chnkr2.r = chnkr2.r + [3;0]; % right rectangle

chnkrs = [chnkr1,chnkr2];
chnkr = merge(chnkrs);

% plot geometry and data

figure(1)
clf
plot(chnkr,'-b')
hold on
quiver(chnkr,'r')
axis equal

%


% solve and visualize the solution

% build CFIE

fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'C',1);

start = tic; sysmat = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = 0.5*eye(chnkr.k*chnkr.nch) + sysmat;

% get the boundary data for a source located at the point above

kerns = @(s,t) chnk.helm2d.kern(zk,s,t,'s');
targs = chnkr.r; targs = reshape(targs,2,chnkr.k*chnkr.nch);
targstau = tangents(chnkr); 
targstau = reshape(targstau,2,chnkr.k*chnkr.nch);

srcinfo = []; srcinfo.r = src0; targinfo = []; targinfo.r = targs;
kernmats = kerns(srcinfo,targinfo);
ubdry = -kernmats*strengths;

rhs = ubdry; rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1);
yl = rmax(2)-rmin(2);
nplot = 300;
xtarg = linspace(-6,6,nplot);
ytarg = linspace(-6,6,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

%

start = tic; in = chunkerinterior(chnkr,targets); t1 = toc(start);
out = ~in;

fprintf('%5.2e s : time to find points in domain\n',t1)

% compute layer potential at points for visualization

optseval = []; optseval.eps = 1e-3;
start = tic;
uscat = chunkerkerneval(chnkr,fkern,sol,targets(:,out),optseval); 
t1 = toc(start);

fprintf('%5.2e s : time to evaluate kernel\n',t1)

srcinfo = []; srcinfo.r = src0; targinfo = []; targinfo.r = targets(:,out);
uin = kerns(srcinfo,targinfo)*strengths;
utot = uscat(:)+uin(:);

%

maxin = max(abs(uin(:)));
maxsc = max(abs(uscat(:)));
maxtot = max(abs(utot(:)));

maxu = max(max(maxin,maxsc),maxtot);

%

figure(2)
clf
subplot(1,3,1)
zztarg = nan(size(xxtarg));
zztarg(out) = uin;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'g')
axis equal
axis tight
colormap(redblue)
caxis([-maxu/5,maxu/5])
title('$u_{in}$','Interpreter','latex','FontSize',24)


subplot(1,3,2)
zztarg = nan(size(xxtarg));
zztarg(out) = uscat;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'g')
axis equal
axis tight
colormap(redblue)
caxis([-maxu/5,maxu/5])
title('$u_{scat}$','Interpreter','latex','FontSize',24)

subplot(1,3,3)
zztarg = nan(size(xxtarg));
zztarg(out) = utot;
h=pcolor(xxtarg,yytarg,imag(zztarg));
set(h,'EdgeColor','none')
hold on
plot(chnkr,'g')
axis equal
axis tight
colormap(redblue)
caxis([-maxu/5,maxu/5])
title('$u_{tot}$','Interpreter','latex','FontSize',24)

