%TEST_CHUNKERIN tests the routines for checking whether a point is 
% inside a domain or not
% 

seed = 8675309;
rng(seed);
addpaths_loc();

doadap = false;

% geometry parameters and construction


cparams = [];
cparams.eps = 1.0e-4;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.5;
start = tic; chnkr = chunkfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

nxdir = 20;
rmin = min(chnkr);
rmax = max(chnkr);
xgrid = linspace(rmin(1),rmax(1),nxdir);
ygrid = linspace(rmin(2),rmax(2),nxdir);
[xx,yy] = meshgrid(xgrid,ygrid);
nt = length(xx(:)); targs = zeros(2,nt); 
targs(1,:) = xx(:); targs(2,:) = yy(:);

opts = [];
opts.verb = false; opts.quadgkparams = {'AbsTol',1e-7,'RelTol',1e-7};
opts.gausseps = 1e-6;
start = tic; in = chunkerin(chnkr,targs,opts); toc(start)

%%

figure(1)
clf
hold off
scatter(targs(1,in),targs(2,in),'go')
hold on
scatter(targs(1,~in),targs(2,~in),'rx')
plot(chnkr,'b','LineWidth',2)
axis equal
