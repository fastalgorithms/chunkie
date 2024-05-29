%CHUNKERMATAPPLYTEST test the routines for a matrix free apply of the system
%matrix
%
% 

clearvars; clear all;
addpaths_loc();

seed = 8675309;
rng(seed);

% chunker geometry parameters and construction

cparams = [];
cparams.eps = 1.0e-4;
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.5;
chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams);

% scalar chunker test
fkern = kernel('lap','d');

start = tic; sysmat = chunkermat(chnkr,fkern); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = -0.5*eye(chnkr.k*chnkr.nch) + sysmat;

fkernsrc = kernel('lap','s');
sources = [1;1];
strengths = [1];
% eval u on bdry

srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = chnkr.r(:,:); 
targinfo.d = chnkr.d(:,:);
dens = fkernsrc.fmm(1e-12,srcinfo,targinfo,strengths);

udense = sys*dens;
start = tic; cormat = chunkermat(chnkr,fkern,struct("corrections",true)); 
toc(start)
sysapply = @(sigma) -0.5*sigma + chunkermatapply(chnkr,fkern,sigma,cormat);
start = tic; u = sysapply(dens); t1 = toc(start);

fprintf('%5.2e s : time for matrix free apply\n',t1)
relerr = norm(udense-u)/norm(udense);
fprintf('relative apply error %5.2e\n',relerr);
assert(relerr < 1e-13)

start = tic; sol1 = sys\dens; toc(start)

start = tic; sol2 = gmres(sysapply, dens, [], 1e-14, 100); toc(start)
relerr = norm(sol1-sol2)/norm(sol1);
fprintf('relative solve error %5.2e\n',relerr);
assert(relerr < 1e-13)

% vector valued chunker test
nregions = 2;
ks = [1.1;2.1]*30;
coefs = [1.0;1.0];
ncurve = 1;
cs(1,1:ncurve) = 1;
cs(2,1:ncurve) = 2;
opts = [];
opts.bdry_data_type = 'point sources';

sources = cell(1,2);
sources{1} = [0.0;0.1];
sources{2} = [3.0;-3.2];
charges{1} = (1.2+1j)*10;
charges{2} = (1+0j)*10;

opts.sources = cell(1,2);
opts.sources{1} = sources{1};
opts.sources{2} = sources{2};
opts.charges = cell(1,2);
opts.charges{1} = charges{1};
opts.charges{2} = charges{2};
[kerns,bdry_data] = chnk.helm2d.transmission_helper(chnkr,ks,cs,coefs,opts);

start = tic; sysmat = chunkermat(chnkr,kerns);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = eye(chnkr.k*chnkr.nch*2) + sysmat;

udense = sys*bdry_data;
cormat = chunkermat(chnkr,kerns,struct("corrections",true));
sysapply = @(sigma) sigma + chunkermatapply(chnkr,kerns,sigma,cormat);

start = tic; u = sysapply(bdry_data); t1 = toc(start);
fprintf('%5.2e s : time for matrix free apply\n',t1)

relerr = norm(udense-u)/norm(udense);
fprintf('relative apply error %5.2e\n',relerr);
assert(relerr < 1e-13)

sol1 = sys\bdry_data;
sol2 = gmres(sysapply, bdry_data, [], 1e-12, 1000);

relerr = norm(sol1-sol2)/norm(sol1);
fprintf('relative solve error %5.2e\n',relerr);
assert(relerr < 1e-10)


% setup chunkgraph 

nverts = 3; 
verts = exp(1i*2*pi*(0:(nverts-1))/nverts);
verts = [real(verts);imag(verts)];

iind = 1:nverts;
jind = 1:nverts;

iind = [iind iind];
jind = [jind jind + 1];
jind(jind>nverts) = 1;
svals = [-ones(1,nverts) ones(1,nverts)];
edge2verts = sparse(iind,jind,svals,nverts,nverts);

amp = 0.1;
frq = 2;
fchnks    = cell(1,size(edge2verts,1));
for icurve = 1:size(edge2verts,1)
    fchnks{icurve} = @(t) sinearc(t,amp,frq);
end
cparams = [];
cparams.nover = 2;
[cgrph] = chunkgraph(verts,edge2verts,fchnks,cparams);

vstruc = procverts(cgrph);
rgns = findregions(cgrph);
cgrph = balance(cgrph);

% scalar chunkgraph test
fkern = -2*kernel('lap','d');

start = tic; sysmat = chunkermat(cgrph,fkern); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = eye(cgrph.npt) + sysmat;

fkernsrc = kernel('lap','s');
sources = [1;1];
strengths = [1];
srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = cgrph.r(:,:); 
targinfo.d = cgrph.d(:,:);
dens = fkernsrc.fmm(1e-12,srcinfo,targinfo,strengths);

udense = sys*dens;
cormat = chunkermat(cgrph,fkern,struct("corrections",true));
sysapply = @(sigma) sigma + chunkermatapply(cgrph,fkern,sigma,cormat);

start = tic; u = sysapply(dens); t1 = toc(start);
fprintf('%5.2e s : time for matrix free apply\n',t1)

relerr = norm(udense-u)/norm(udense);
fprintf('relative apply error %5.2e\n',relerr);
assert(relerr < 1e-13)

sol1 = sys\dens;

sol2 = gmres(sysapply, dens, [], 1e-14, 100);
relerr = norm(sol1-sol2)/norm(sol1);
fprintf('relative solve error %5.2e\n',relerr);
assert(relerr < 1e-13)

% vectorvalued chunkgraph test
nregions = 2;
ks = [1.1;2.1]*30;
coefs = [1.0;1.0];
ncurve = 3;
cs(1,1:ncurve) = 1;
cs(2,1:ncurve) = 2;
opts = [];
opts.bdry_data_type = 'point sources';

sources = cell(1,2);
sources{1} = [0.0;0.1];
sources{2} = [3.0;-3.2];
charges{1} = (1.2+1j)*10;
charges{2} = (1+0j)*10;

opts.sources = cell(1,2);
opts.sources{1} = sources{1};
opts.sources{2} = sources{2};
opts.charges = cell(1,2);
opts.charges{1} = charges{1};
opts.charges{2} = charges{2};
[kerns,bdry_data] = chnk.helm2d.transmission_helper(cgrph,ks,cs,coefs,opts);

start = tic; sysmat = chunkermat(cgrph,kerns); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = eye(size(sysmat,1)) + sysmat;


udense = sys*bdry_data;
cormat = chunkermat(cgrph,kerns,struct("corrections",true));
sysapply = @(sigma) sigma + chunkermatapply(cgrph,kerns,sigma,cormat);
start = tic; u = sysapply(bdry_data); t1 = toc(start);
fprintf('%5.2e s : time for matrix free apply\n',t1)

relerr = norm(udense-u)/norm(udense);
fprintf('relative apply error %5.2e\n',relerr);
assert(relerr < 1e-13)

sol1 = sys\bdry_data;
sol2 = gmres(sysapply, bdry_data, [], 1e-12, 200);

relerr = norm(sol1-sol2)/norm(sol1);
fprintf('relative solve error %5.2e\n',relerr);
assert(relerr < 1e-10)


function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
