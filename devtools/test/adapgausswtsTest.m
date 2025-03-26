
%ADAPGAUSSWTSTEST
%
% define geometry and test adaptive integration routine

clearvars; close all;
iseed = 8675309;
rng(iseed);

zk = randn() + 1i*randn();

addpaths_loc();

cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 2;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.25;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

fprintf('%5.2e s : time to build geo\n',t1)

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

% targets

nt = 3;
ts = 0.0+2*pi*rand(nt,1);
targets = starfish(ts,narms,amp);
targets = targets.*repmat(rand(1,nt),2,1);

% plot geo and sources

% xs = chnkr.r(1,:,:); xmin = min(xs(:)); xmax = max(xs(:));
% ys = chnkr.r(2,:,:); ymin = min(ys(:)); ymax = max(ys(:));
% 
% figure(1)
% clf
% hold off
% plot(chnkr)
% hold on
% scatter(sources(1,:),sources(2,:),'o')
% scatter(targets(1,:),targets(2,:),'x')
% axis equal 

%

% build layer potential matrix with GGQ routine for comparison

%fkern = @(s,t) chnk.lap2d.kern(s,t,'S');
fkern = @(s,t) chnk.helm2d.kern(zk,s,t,'D');
start = tic; mat1 = chunkermat(chnkr,fkern);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

%

% use adaptive routine on a chunk to generate new entries for 
% a neighbor interaction

% normally, this would be called inside a matrix/sparse matrix builder
% we need to grab a lot of info...

nch = chnkr.nch;
j = 1;
ich = chnkr.adj(1,j);

opdims = zeros(2,1);
opdims(1) = 1; opdims(2) = 1;

r = chnkr.r; d = chnkr.d; d2 = chnkr.d2; n =chnkr.n;
k = chnkr.k; [t,w] = lege.exps(k);
bw = lege.barywts(k);

k2 = max(27,k+1);
[t2,w2] = lege.exps(k2);

rt = r(:,:,ich);
dt = d(:,:,ich);
d2t = d2(:,:,ich);
nt = n(:,:,ich);
dtlen = sqrt(sum(dt.^2,1));
taut = dt./dtlen;

opts = [];
opts.eps = 1e-5;

ntimes = 100;
start = tic;
for i = 1:ntimes
    [mat,maxrecs,numints,iers] = chnk.adapgausswts(r,d,n,d2,[],t,bw,j, ...
        rt,dt,nt,d2t,[],fkern,opdims,t2,w2,opts);
end
t1 = toc(start);

fprintf('speed %5.2e\n',ntimes/t1);


% compare

istart = opdims(1)*(ich-1)*k + 1;
iend = opdims(1)*ich*k;
jstart = opdims(1)*(j-1)*k + 1;
jend = opdims(1)*j*k;

matcomp = mat1(istart:iend,jstart:jend);

assert(norm(matcomp-mat,'inf') < 1e-11);
