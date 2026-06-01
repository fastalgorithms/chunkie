
chunkermat_axissymlap_DirichletTest0();

function chunkermat_axissymlap_DirichletTest0()

iseed = 8675309;
rng(iseed);

%addpaths_loc();

type = 'chnkr-star';
type = 'chnkr-torus';
type = 'cgrph';
%type = 'cgrph-sphere';

pref = [];
pref.k = 16;
ns = 100;
nt = 200;
%ns = 1;
%nt = 1;

maxchunklen = 0.5;

[chnkr, targets, sources] = get_geometry(type, pref, nt, ns, maxchunklen);

%targets = [0.4,0.3;-0.2,0.2];
%sources = [0.01;1.4];

fprintf('Done building geometry\n');

% source strengths
strengths = randn(ns, 1);

% targets

% Plot everything

if (false)
figure(1)
clf
hold off
plot(chnkr, 'k.');
hold on
quiver(chnkr);
scatter(sources(1,:), sources(2,:), 'o','green','filled')
scatter(targets(1,:), targets(2,:), 'x','blue')
axis equal

end

% Set up kernels
D      = kernel('axissymlaplace','d');
S      = kernel('axissymlaplace','s');

K = -1/(2*pi^2)*D;
Keval = K;

opts = [];
opts.rcip = true;
opts.nsub_or_tol = 40;
npts = chnkr.npt;
nsys = K.opdims(1)*npts;
tic, A = chunkermat(chnkr, K, opts) + eye(nsys); toc

srcinfo  = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.n = chnkr.n(:,:);
kernmats = S.eval(srcinfo, targinfo);
ubdry = kernmats*strengths;
rhs = ubdry;

start = tic;
sol = A\rhs;
%sol = gmres(A, rhs, [], 1e-14, 200);
t1 = toc(start);

% Compute exact solution
srcinfo  = []; srcinfo.r  = sources;
targinfo = []; targinfo.r = targets;
kernmatstarg = S.eval(srcinfo, targinfo);
utarg = kernmatstarg*strengths;


% Compute solution using chunkerkerneval
% evaluate at targets and compare

opts.forcesmooth = false;
opts.verb = false;
opts.quadkgparams = {'RelTol', 1e-15, 'AbsTol', 1.0e-15};

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

assert(relerr < 1e-10)
assert(relerr2 < 1e-10)

end

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
    sources = 0.2*[cos(ts)';sin(ts)'];
    
    ts = -pi/2 + pi*rand(nt,1);
    targets = 3.0*[cos(ts)'; sin(ts)'];



elseif strcmpi(type,'cgrph-sphere')
    
    
    nverts = 2; 
    verts = exp(-1i*pi/2 + 1i*pi*(0:(nverts-1))/(nverts-1));
    verts = [real(verts);imag(verts)];


    iind = 1:(nverts-1);
    jind = 1:(nverts-1);

    iind = [iind iind];
    jind = [jind jind + 1];
    jind(jind>nverts) = 1;
    svals = [-ones(1,nverts-1) ones(1,nverts-1)];
    edge2verts = sparse(iind, jind, svals, nverts-1, nverts);

    cparams = [];
    cparams.eps = 1.0e-10;
    cparams.nover = 1;
    cparams.ifclosed = false;
    cparams.ta = 0;
    cparams.tb = 1;
    cparams.maxchunklen = maxchunklen;
    narms = 0;
    amp = 0.0;
    
    fchnks{1} = @(t) circle(t);

    
    chnkobj = chunkgraph(verts, edge2verts, fchnks, cparams, pref);
    chnkobj = balance(chnkobj);
    
    ts = -pi/2 + pi*rand(ns, 1);
    sources = starfish(ts, narms, amp);
    sources = 1.2*sources;

    ts = -pi/2 + pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (0.8*repmat(rand(1, nt), 2, 1));   

    
    
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


function [r,d,d2] = circle(t)
xs = cos(-pi/2 + pi*t);
ys = sin(-pi/2 + pi*t);
xp = -pi*ys;
yp = pi*xs;
xpp = -pi*pi*xs;
ypp = -pi*pi*ys;
r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end



