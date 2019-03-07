%TEST_GAUSSID test the routines for integrating over chunks against the
% Green's id for Laplace
%
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

% scalar function on boundary

kernd = @(s,t,sn,tn) glapkern(s,t,sn,tn,'d');
kerns = @(s,t,sn,tn) glapkern(s,t,sn,tn,'s');

dens1 = ones(chunker.k,chunker.nch);

ndims = [1 1];

nxdir = 100;
xs = chunker.chunks(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chunker.chunks(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));
xgrid = linspace(xmin,xmax,nxdir);
ygrid = linspace(ymin,ymax,nxdir);
[xx,yy] = meshgrid(xgrid,ygrid);
nt = length(xx(:)); targs = zeros(2,nt); 
targs(1,:) = xx(:); targs(2,:) = yy(:);

opts.usesmooth=true;
start=tic; d1 = chunkerintkern(chunker,kernd,ndims,dens1,targs,opts); 
toc(start)

if doadap
    opts.usesmooth=false;
    opts.verb=true;
    opts.quadkgparams = {'RelTol',1.0e-2,'AbsTol',1.0e-2};
    start=tic; d12 = chunkerintkern(chunker,kernd,ndims,dens1,targs,opts); 
    toc(start)
    dd2 = reshape(d12,size(xx));
end

dd = reshape(d1,size(xx));

%%

hold off

h = pcolor(xx,yy,1.0*or((abs(dd+1)<1.0e-8),abs(dd)<1.0e-8)); set(h,'EdgeColor','none');
hold on
scatter(xs(:),ys(:))

colorbar

