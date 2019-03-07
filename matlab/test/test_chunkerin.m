%TEST_CHUNKIN tests the routines for checking whether a point is 
% inside a domain or not
% 

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../../mwrap')

doadap = false;

% geometry parameters and construction

cparams.eps = 1.0e-6;
cparams.nchmax = 100000;
cparams.nover = 1;
narms = 5;
amp = 0.5;
chunker = chunkfunc(@(t) starfish(t,narms,amp),cparams);

nxdir = 20;
xs = chunker.chunks(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chunker.chunks(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));
xgrid = linspace(xmin,xmax,nxdir);
ygrid = linspace(ymin,ymax,nxdir);
[xx,yy] = meshgrid(xgrid,ygrid);
nt = length(xx(:)); targs = zeros(2,nt); 
targs(1,:) = xx(:); targs(2,:) = yy(:);

opts = [];
opts.verb = true; opts.quadgkparams = {'AbsTol',1e-3,'RelTol',1e-3};
opts.gausseps = 1e-6;
start = tic; in = chunkerin(chunker,targs,opts); toc(start)

%

hold off
scatter(targs(1,in),targs(2,in),'go')
hold on
scatter(targs(1,~in),targs(2,~in),'rx')
scatter(xs(:),ys(:),'bo')

colorbar

