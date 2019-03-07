%TEST_GAUSSID test the routines for integrating over chunks against the
% Gauss' identity for the double layer potential
%
% 

seed = 8675309;
rng(seed);
addpath('../src')
addpath('../../mwrap')

doadap = true;

% geometry parameters and construction

cparams.eps = 1.0e-3;
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

nxdir = 40;
xs = chunker.chunks(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
ys = chunker.chunks(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));
xgrid = linspace(xmin,xmax,nxdir);
ygrid = linspace(ymin,ymax,nxdir);
[xx,yy] = meshgrid(xgrid,ygrid);
nt = length(xx(:)); targs = zeros(2,nt); 
targs(1,:) = xx(:); targs(2,:) = yy(:);

fprintf('computing Gauss I.D. with smooth rule...\n');
opts.usesmooth=true;
start=tic; d1 = chunkerintkern(chunker,kernd,ndims,dens1,targs,opts); 
toc(start)

if doadap
    fprintf( ...
      'computing Gauss I.D. with adaptive quadrature (may be slow)...\n');
    opts.usesmooth=false;
    opts.verb=false;
    opts.quadkgparams = {'RelTol',1.0e-7,'AbsTol',1.0e-7};
    start=tic; d12 = chunkerintkern(chunker,kernd,ndims,dens1,targs,opts); 
    toc(start)
    dd2 = reshape(d12,size(xx));
end

dd = reshape(d1,size(xx));

%%

figure(1)

hold off

h = pcolor(xx,yy,1.0*or((abs(dd+1)<1.0e-6),abs(dd)<1.0e-6)); set(h,'EdgeColor','none');
hold on
scatter(xs(:),ys(:))

caxis([0 1]);

colorbar

if doadap
    figure(2)

    hold off

    h = pcolor(xx,yy,1.0*or((abs(dd2+1)<1.0e-8),abs(dd2)<1.0e-8)); set(h,'EdgeColor','none');
    hold on
    scatter(xs(:),ys(:))

    caxis([0 1])
    
    colorbar
end

