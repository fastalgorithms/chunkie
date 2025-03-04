%CHUNKERINTERIORTEST tests the routines for checking whether a point is 
% inside a domain or not
% 

seed = 8675309;
rng(seed);

doadap = false;

% geometry parameters and construction


cparams = [];
cparams.eps = 1.0e-9;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.5;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

% nxdir = 20;
% rmin = min(chnkr);
% rmax = max(chnkr);
% xgrid = linspace(rmin(1),rmax(1),nxdir);
% ygrid = linspace(rmin(2),rmax(2),nxdir);
% [xx,yy] = meshgrid(xgrid,ygrid);
% nt = length(xx(:)); targs = zeros(2,nt); 
% targs(1,:) = xx(:); targs(2,:) = yy(:);

nt = 10000;
scal = 2*rand(1,nt);
tr = 2*pi*rand(1,nt);

targs = bsxfun(@times,starfish(tr,narms,amp),scal);


opts = [];
opts.flam = false;
opts.fmm = false;
start = tic; in = chunkerinterior(chnkr,targs,opts); t1 = toc(start);

opts = [];
opts.fmm = false;
opts.flam = true;
start = tic; in2 = chunkerinterior(chnkr,targs,opts); t2 = toc(start);

opts = [];
opts.fmm = true;
opts.flam = false;
start = tic; in3 = chunkerinterior(chnkr,targs,opts); t3 = toc(start);


fprintf('%5.2e s : time for chunkerinterior (no flam)\n',t1);
fprintf('%5.2e s : time for chunkerinterior (with flam)\n',t2);
fprintf('%5.2e s : time for chunkerinterior (with fmm)\n',t3);

assert(all(in(:) == (scal(:) < 1)));
assert(all(in2(:) == (scal(:) < 1)));
assert(all(in3(:) == (scal(:) < 1)));

x1 = linspace(-2,2,100);
[xx,yy] = meshgrid(x1,x1);

opts = [];
opts.fmm = true;
opts.flam = false;
start = tic; in4 = chunkerinterior(chnkr,{x1,x1},opts); t4 = toc(start);
fprintf('%5.2e s : time for chunkerinterior (meshgrid, with fmm)\n',t4);

% Test targets specified as chunker
narms = 3;
amp = 0.1;
chnkr2 = chunkerfunc(@(t) 0.3*starfish(t,narms,amp),cparams,pref); 

opts = [];
opts.fmm = true;
opts.flam = false;
start = tic; in5 = chunkerinterior(chnkr,chnkr2,opts); t5 = toc(start);
fprintf('%5.2e s : time for chunkerinterior (chunker, with fmm)\n',t5);
assert(all(in5(:) == 1));

% test axissym option 
chnkr = chunkerfunc(@(t) starfish(t),struct('ta',-pi/2,'tb',pi/2,'ifclosed',0));
nt = 1000;
ttarg = -pi/2+pi*rand(nt,1); scal = 2*rand(1,nt);
targs = starfish(ttarg).*scal;

in = chunkerinterior(chnkr,targs,struct('axissym',true));
assert(all(in(:) == (scal(:) <= 1)));

%% 

% a stress test. previously this triggered a bug for the old version which
% used the normals to determine interior points

amp = 0.25;
scale = .3;

ctr = [-2;-1.6];
chnkr_int= chunkerfunc(@(t) starfish(t,3,amp,ctr,pi/4,scale)); 
chnkr_int = sort(reverse(chnkr_int));

a = max(vecnorm(chnkr_int.r(:,:)))*1.01;
chnkr_ext = chunkerfunc(@(t) [a*cos(t(:).'); a*sin(t(:).')]); 
chnkr = merge([chnkr_ext,chnkr_int]);

% return

L = max(abs(chnkr.r),[],"all");
x1 = linspace(-L,L,100);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).']; 
xv=[chnkr_ext.r(1,:),nan,chnkr_int.r(1,:)];
yv=[chnkr_ext.r(2,:),nan,chnkr_int.r(2,:)];
%tic; in0 = inpolygon(xx(:),yy(:),xv,yv); toc
tt = atan2(yy-ctr(2),xx-ctr(1)) + 2*pi;
st = starfish(tt(:),3,amp,[0;0],pi/4,scale);
ss2 = reshape(st(1,:).^2 + st(2,:).^2,size(xx));
in0 = and((xx.^2+yy.^2)<a^2, ...
    (xx-ctr(1)).^2 + (yy-ctr(2)).^2 > ss2);
in0 = in0(:);
tic; in = chunkerinterior(chnkr,{x1,x1}); toc

% clf
% plot(chnkr,'g-o')
% hold on
% scatter(xx(in(:)),yy(in(:)),'bo')
% scatter(xx(~in(:)),yy(~in(:)),'rx')

assert(all(in0 == in));

