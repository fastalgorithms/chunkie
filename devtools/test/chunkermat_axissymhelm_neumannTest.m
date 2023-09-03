clearvars; close all;
iseed = 8675309;
rng(iseed);

addpaths_loc();

zk = 1j*10.1;

type = 'chnkr-star';
% type = 'chnkr-torus';

irep = 'rpcomb';
irep = 'sk';

pref = [];
pref.k = 16;
ns = 10;
nt = 100;
ppw = 80;   % points per wavelength;
maxchunklen = pref.k/ppw/real(zk)*2*pi;
maxchunklen = 0.5;

[chnkr, sources, targets] = get_geometry(type, pref, ns, nt, maxchunklen);
wts = chnkr.wts; wts = wts(:);

l2scale = false;
fprintf('Done building geometry\n');

% source strengths
strengths = randn(ns, 1);

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


% For solving exterior the Neumann boundary value problem, we use the
% following integral equation
%
% u = \beta(S_{k} + i \alpha D_{k} S_{ik}) \sigma
%
% with \beta = -1.0/(0.5 + 1i*0.25*alpha)
% 
% On imposing the boundary conditions, we get the following integral 
% equation
%
%  du/dn = I + \beta S_{k}'[\sigma] + 
%          i\beta \alpha(D'_{k} - D'_{ik})(S_{ik}[\sigma]) + 
%          i\beta \alpha(S_{ik}')^2 [\sigma];
% 
% Setting -S_{ik}[\sigma] + \tau = 0, and - S_{ik}'[\sigma] + \mu=0, 
% we have the following system of integral equations

% Set up kernels
Skp    = kernel('axissymhelm', 'sprime', zk);
Sk     = kernel('axissymhelm', 's', zk);
Dk     = kernel('axissymhelm', 'd', zk);

Z = kernel.zeros();


if strcmpi(irep,'rpcomb')
    Sik    = kernel('axissymhelm', 's', 1i*zk);
    Sikp   = kernel('axissymhelm', 'sprime', 1i*zk);
    Dkdiff = kernel('axissymhelmdiff', 'dprime', [zk 1i*zk]);
    alpha = 1;
    c1 = -1/(0.5 + 1i*alpha*0.25);
    c2 = -1i*alpha/(0.5 + 1i*alpha*0.25);
    c3 = -1;
    K = [ c1*Skp  c2*Dkdiff c2*Sikp ;
       c3*Sik  Z        Z        ;
       c3*Sikp Z        Z        ];
    K = kernel(K);
    Keval = c1*kernel([Sk 1i*alpha*Dk Z]);
else
    K = -2*Skp;
    Keval = -2*Sk;
end

% Set up boundary data

srcinfo  = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.n = chnkr.n(:,:);
kernmats = Skp.eval(srcinfo, targinfo);
ubdry = kernmats*strengths;

npts = chnkr.npt;
nsys = K.opdims(1)*npts;
rhs = zeros(nsys, 1);


if(l2scale)
    rhs(1:K.opdims(1):end) = ubdry.*sqrt(wts);
else
    rhs(1:K.opdims(1):end) = ubdry;
end

% Form matrix
opts = [];
opts.l2scale = l2scale;
tic, A = chunkermat(chnkr, K, opts) + eye(nsys); toc
start = tic;
sol = gmres(A, rhs, [], 1e-14, 100);
t1 = toc(start);

% Compute exact solution
srcinfo  = []; srcinfo.r  = sources;
targinfo = []; targinfo.r = targets;
kernmatstarg = Sk.eval(srcinfo, targinfo);
utarg = kernmatstarg*strengths;

% Compute solution using chunkerkerneval
% evaluate at targets and compare

opts.usesmooth = false;
opts.verb = false;
opts.quadkgparams = {'RelTol', 1e-16, 'AbsTol', 1.0e-16};

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
Dsol = chunkerkerneval(chnkrtotal, Keval, sol, targets, opts);
t2 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n', t2)


wchnkr = chnkrtotal.wts;
wchnkr = repmat(wchnkr(:).', K.opdims(1), 1);
relerr  = norm(utarg-Dsol) / (sqrt(chnkrtotal.nch)*norm(utarg));
relerr2 = norm(utarg-Dsol, 'inf') / dot(abs(sol(:)), wchnkr(:));
fprintf('relative frobenius error %5.2e\n', relerr);
fprintf('relative l_inf/l_1 error %5.2e\n', relerr2);

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

    chnkobj = chunkgraph(verts, edge2verts, fchnks, cparams, pref);
    chnkobj = balance(chnkobj);
    
    ts = 0.0+2*pi*rand(ns,1);
    sources = 3.0*[cos(ts)';sin(ts)'];
    
    ts = 0.0+2*pi*rand(nt,1);
    targets = 0.2*[cos(ts)'; sin(ts)'];


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
    sources = 0.5*sources;

    ts = -pi/2 + pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (1 + 0.5*repmat(rand(1, nt), 2, 1));   

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
    sources = 0.5*sources;
    sources(1,:) = sources(1,:) + ctr(1);

    ts = -pi/2 + pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (1 + 0.5*repmat(rand(1, nt), 2, 1));    
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

