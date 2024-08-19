clearvars; close all;
iseed = 8675309;
rng(iseed);

addpaths_loc();

zk = 10.1;

type = 'cgrph';
% type = 'chnkr-torus';

pref = [];
pref.k = 16;
ns = 10;
nt = 10;
ppw = 80;   % points per wavelength;
maxchunklen = pref.k/ppw/real(zk)*2*pi;
maxchunklen = 0.5;

[chnkr, sources, targets] = get_geometry(type, pref, ns, nt, maxchunklen);
wts = chnkr.wts; wts = wts(:);

l2scale = false;
fprintf('Done building geometry\n');

% source strengths
strengths_in = randn(ns, 1);
strengths_out = randn(nt, 1);


% targets


% Plot everything

figure(1)
clf
hold off
plot(chnkr)
hold on
scatter(sources(1,:), sources(2,:), 'o')
scatter(targets(1,:), targets(2,:), 'x')
axis equal


% For solving the transmission boundary value problem, we
% use the repesentation 
% u_{i} = D_{k_{i}} [\sigma] - S_{k_{i}}[\tau]
% 
% and impose jumps in u and du/dn as exterior-interior respectively
%

% Set up kernels
Skdiff     = kernel('axissymhelmdiff', 's', [zk 1i*zk]);
Dkdiff     = kernel('axissymhelmdiff', 'd', [zk 1i*zk]);
Skpdiff    = kernel('axissymhelmdiff', 'sprime', [zk 1i*zk]);
Dkpdiff    = kernel('axissymhelmdiff', 'dprime', [zk 1i*zk]);

K = [Dkdiff  -1*Skdiff;
     Dkpdiff    -1*Skpdiff];
K = kernel(K);
Dk    = kernel('axissymhelm', 'd', zk);
Sk    = kernel('axissymhelm', 's', zk);
Skp    = kernel('axissymhelm', 'sprime', zk);

Dik    = kernel('axissymhelm', 'd', 1i*zk);
Sik    = kernel('axissymhelm', 's', 1i*zk);
Sikp    = kernel('axissymhelm', 'sprime', 1i*zk);

Kouteval =  kernel([Dk, -1*Sk]);
Kineval =  kernel([Dik, -1*Sik]);


% Set up boundary data

srcinfo  = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.n = chnkr.n(:,:);
kernmats = Sk.eval(srcinfo, targinfo);
kernmatsp = Skp.eval(srcinfo, targinfo);
ubdry = kernmats*strengths_in;
dudnbdry = kernmatsp*strengths_in;



srcinfo  = []; srcinfo.r = targets; 
targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.n = chnkr.n(:,:);
kernmats = Sik.eval(srcinfo, targinfo);
kernmatsp = Sikp.eval(srcinfo, targinfo);
ubdry = ubdry - kernmats*strengths_out;
dudnbdry = dudnbdry - kernmatsp*strengths_out;


npts = chnkr.npt;
nsys = K.opdims(1)*npts;
rhs = zeros(nsys, 1);


if(l2scale)
    rhs(1:2:end) = ubdry.*sqrt(wts);
    rhs(2:2:end) = dudnbdry.*sqrt(wts);
else
    rhs(1:2:end) = ubdry;
    rhs(2:2:end) = dudnbdry;
end

% Form matrix
opts = [];
opts.l2scale = l2scale;
tic, A = chunkermat(chnkr, K, opts) + eye(nsys); toc
start = tic;
sol = gmres(A, rhs, [], 1e-14, 200);
t1 = toc(start);

% Compute exact solution
srcinfo  = []; srcinfo.r  = sources;
targinfo = []; targinfo.r = targets;
kernmatstarg = Sk.eval(srcinfo, targinfo);
utarg_out = kernmatstarg*strengths_in;

% Compute exact solution
srcinfo  = []; srcinfo.r  = targets;
targinfo = []; targinfo.r = sources;
kernmatstarg = Sik.eval(srcinfo, targinfo);
utarg_in = kernmatstarg*strengths_out;



% Compute solution using chunkerkerneval
% evaluate at targets and compare

opts.forceadap = true;
opts.verb = false;
opts.quadkgparams = {'RelTol', 1e-8, 'AbsTol', 1.0e-8};


if(l2scale)
    wts_rep = repmat(wts(:).', K.opdims(1),1);
    wts_rep = wts_rep(:);
    sol = sol./sqrt(wts_rep);
end

if isa(chnkr, 'chunkgraph')
    % Collapse cgrph into chnkrtotal
    chnkrs = chnkr.echnks;
    chnkrtotal = merge(chnkrs);
else
    chnkrtotal = chnkr;
end



start = tic;
Dsol_out = chunkerkerneval(chnkrtotal, Kouteval, sol, targets, opts);
Dsol_in =  chunkerkerneval(chnkrtotal, Kineval, sol, sources, opts);
t2 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n', t2)


wchnkr = chnkrtotal.wts;
wchnkr = repmat(wchnkr(:).', K.opdims(1), 1);
relerr  = norm(utarg_out-Dsol_out) / (sqrt(chnkrtotal.nch)*norm(utarg_out));
relerr2 = norm(utarg_out-Dsol_out, 'inf') / dot(abs(sol(:)), wchnkr(:));
fprintf('relative frobenius error in exterior %5.2e\n', relerr);
fprintf('relative l_inf/l_1 error in exterior %5.2e\n', relerr2);

relerr  = norm(utarg_in-Dsol_in) / (sqrt(chnkrtotal.nch)*norm(utarg_in));
relerr2 = norm(utarg_in-Dsol_in, 'inf') / dot(abs(sol(:)), wchnkr(:));
fprintf('relative frobenius error in interior %5.2e\n', relerr);
fprintf('relative l_inf/l_1 error in interior %5.2e\n', relerr2);


return




% Test fast direct solver interfaces

% build sparse tridiag part


opts.nonsmoothonly = true;
opts.rcip = true;
start = tic; spmat = chunkermat(chnkr, K, opts); t1 = toc(start);
fprintf('%5.2e s : time to build tridiag\n',t1)


spmat = spmat + speye(nsys);

% test matrix entry evaluator
start = tic; 
opdims = K.opdims;
sys2 = chnk.flam.kernbyindex(1:nsys, 1:nsys, chnkr, K, opdims, ...
    spmat, l2scale);


t1 = toc(start);

fprintf('%5.2e s : time for mat entry eval on whole mat\n',t1)

err2 = norm(sys2-A,'fro')/norm(A,'fro');
fprintf('%5.2e   : fro error of build \n',err2);

% test fast direct solver
opts.ifproxy = false;
F = chunkerflam(chnkr,K,1.0,opts);

start = tic; sol2 = rskelf_sv(F,rhs); t1 = toc(start);

if(l2scale)
    wts_rep = repmat(wts(:).', K.opdims(1),1);
    wts_rep = wts_rep(:);
    sol2 = sol2./sqrt(wts_rep);
end


fprintf('%5.2e s : time for rskelf_sv \n',t1)

err = norm(sol-sol2,'fro')/norm(sol,'fro');

fprintf('difference between fast-direct and iterative %5.2e\n',err)


function [chnkobj, sources, targets] = get_geometry(type, pref, ns, nt, maxchunklen)

if nargin == 0
    type = 'chnkr';
end

if nargin <= 1
    pref = [];
    pref.k = 16;
end

if nargin <= 2
    ns = 10;
end

if nargin <= 3
    nt = 10;
end


if nargin <= 4
    maxchunklen = 1.0;
end
    

if strcmpi(type, 'cgrph')
    
    
    nverts = 3; 
    verts = exp(-1i*pi/2 + 1i*pi*(0:(nverts-1))/(nverts-1));
    verts = [real(verts);imag(verts)];


    iind = 1:(nverts-1);
    jind = 1:(nverts-1);

    iind = [iind iind];
    jind = [jind jind + 1];
    jind(jind>nverts) = 1;
    svals = [-ones(1,nverts-1) ones(1,nverts-1)];
    edge2verts = sparse(iind, jind, svals, nverts-1, nverts);

    amp = 0.1;
    frq = 2;
    fchnks    = cell(1,size(edge2verts,1));
    for icurve = 1:size(edge2verts,1)
        fchnks{icurve} = @(t) sinearc(t, amp, frq);
    end
    cparams = [];
    cparams.nover = 2;
    cparams.maxchunklen = maxchunklen;
    cparams.ta = 0;
    cparams.tb = 1;

    chnkobj = chunkgraph(verts, edge2verts, fchnks, cparams, pref);
    chnkobj = balance(chnkobj);
       
    ts = -pi/2 + pi*rand(ns,1);
    sources = 0.2*[cos(ts)'; sin(ts)'];

    ts = -pi/2 + pi*rand(nt,1);
    targets = 3.0*[cos(ts)';sin(ts)'];


elseif strcmpi(type,'chnkr-star')
    cparams = [];
    cparams.eps = 1.0e-10;
    cparams.nover = 1;
    cparams.ifclosed = false;
    cparams.ta = -pi/2;
    cparams.tb = pi/2;
    cparams.maxchunklen = maxchunklen;
    narms = 0;
    amp = 0.25;
    chnkobj = chunkerfunc(@(t) starfish(t, narms, amp), cparams, pref); 
    chnkobj = sort(chnkobj);
    
    ts = -pi/2 + pi*rand(ns, 1);
    sources = starfish(ts, narms, amp);
    sources = sources .* (0.7 - 0.4.*repmat(rand(1,ns), 2, 1));

    ts = -pi/2 + pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (1.2 + 3.0*repmat(rand(1, nt), 2, 1));   

elseif strcmpi(type,'chnkr-torus')
    cparams = [];
    cparams.eps = 1.0e-10;
    cparams.nover = 1;
    cparams.ifclosed = true;
    cparams.ta = 0;
    cparams.tb = 2*pi;
    cparams.maxchunklen = maxchunklen;
    narms = 0;
    amp = 0.25;
    ctr = [3 0];
    chnkobj = chunkerfunc(@(t) starfish(t, narms, amp, ctr), cparams, pref); 
    chnkobj = sort(chnkobj);
    
    ts = -pi/2 + pi*rand(ns, 1);
    sources = starfish(ts, narms, amp);
    sources = sources .* (0.7 - 0.4.*repmat(rand(1,ns), 2, 1));
    sources(1,:) = sources(1,:) + ctr(1);

    ts = -pi/2 + pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (1.2 + 3.0*repmat(rand(1, nt), 2, 1));    
    targets(1,:) = targets(1,:) + ctr(1);
end



end


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

