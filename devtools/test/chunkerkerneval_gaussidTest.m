%CHUNKERKERNEVAL_GAUSSIDTEST test the routines for integrating over 
% chunks against Gauss' identity for the double layer potential
%
% 

clearvars; close all;
seed = 8675309;
rng(seed);
addpaths_loc();

doadap = true;

% geometry parameters and construction
 
cparams = [];
cparams.eps = 1.0e-4;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.5;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);


% scalar function on boundary

kernd = @(s,t) chnk.lap2d.kern(s,t,'d');
kerns = @(s,t) chnk.lap2d.kern(s,t,'s');

dens1 = ones(chnkr.k,chnkr.nch);

nxdir = 40;

rmin = min(chnkr);
rmax = max(chnkr);
xgrid = linspace(rmin(1),rmax(1),nxdir);
ygrid = linspace(rmin(2),rmax(2),nxdir);
[xx,yy] = meshgrid(xgrid,ygrid);
nt = length(xx(:)); targs = zeros(2,nt); 
targs(1,:) = xx(:); targs(2,:) = yy(:);

fprintf('computing Gauss I.D. with smooth rule...\n');
opts.forcesmooth=true;
opts.flam = true;
start=tic; d1 = chunkerkerneval(chnkr,kernd,dens1,targs,opts); 
toc(start)



if doadap
    fprintf( ...
      'computing Gauss I.D. with adaptive quadrature (may be slow)...\n');
    opts.forcesmooth=false;
    opts.forceadap = false;
    opts.quadkgparams = {'RelTol',1.0e-7,'AbsTol',1.0e-7};
    opts.fac = 1.0;
    opts.flam = true;
    start=tic; d12 = chunkerkerneval(chnkr,kernd,dens1,targs,opts); 
    toc(start)
    assert(all(or(abs(d12(:))<1e-6,abs(d12(:)+1)<1e-6)));
    dd2 = reshape(d12,size(xx));
end

dd = reshape(d1,size(xx));

%

figure(1)
clf
hold off

h = pcolor(xx,yy,log10(min(abs(dd+1),abs(dd)))); set(h,'EdgeColor','none');
hold on
plot(chnkr)

caxis([-16,0])

colorbar

if doadap
    figure(2)
    clf 
    hold off

    h = pcolor(xx,yy,log10(min(abs(dd2+1),abs(dd2)))); set(h,'EdgeColor','none');
    hold on
    plot(chnkr)

    caxis([-16,0])
    
    colorbar
end

