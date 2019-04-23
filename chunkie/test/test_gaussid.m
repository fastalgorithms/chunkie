%TEST_GAUSSID test the routines for integrating over chunks against the
% Gauss' identity for the double layer potential
%
% 

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
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

% scalar function on boundary

kernd = @(s,t,sn,tn) glapkern(s,t,sn,tn,'d');
kerns = @(s,t,sn,tn) glapkern(s,t,sn,tn,'s');

dens1 = ones(chnkr.k,chnkr.nch);

opdims = [1 1];

nxdir = 40;

rmin = min(chnkr);
rmax = max(chnkr);
xgrid = linspace(rmin(1),rmax(1),nxdir);
ygrid = linspace(rmin(2),rmax(2),nxdir);
[xx,yy] = meshgrid(xgrid,ygrid);
nt = length(xx(:)); targs = zeros(2,nt); 
targs(1,:) = xx(:); targs(2,:) = yy(:);

fprintf('computing Gauss I.D. with smooth rule...\n');
opts.usesmooth=true;
start=tic; d1 = chunkerintkern(chnkr,kernd,opdims,dens1,targs,opts); 
toc(start)

if doadap
    fprintf( ...
      'computing Gauss I.D. with adaptive quadrature (may be slow)...\n');
    opts.usesmooth=false;
    opts.verb=false;
    opts.quadkgparams = {'RelTol',1.0e-7,'AbsTol',1.0e-7};
    start=tic; d12 = chunkerintkern(chnkr,kernd,opdims,dens1,targs,opts); 
    toc(start)
    dd2 = reshape(d12,size(xx));
end

dd = reshape(d1,size(xx));

%%

figure(1)
clf
hold off

h = pcolor(xx,yy,1.0*or((abs(dd+1)<1.0e-6),abs(dd)<1.0e-6)); set(h,'EdgeColor','none');
hold on
plot(chnkr)

caxis([0 1]);

colorbar

if doadap
    figure(2)
    clf 
    hold off

    h = pcolor(xx,yy,1.0*or((abs(dd2+1)<1.0e-8),abs(dd2)<1.0e-8)); set(h,'EdgeColor','none');
    hold on
    plot(chnkr)

    caxis([0 1])
    
    colorbar
end

