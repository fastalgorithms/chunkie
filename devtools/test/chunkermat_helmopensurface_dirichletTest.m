% Solve the Helmholtz Dirichlet problem on an open curve, we use the
% following integral representation
%
% u = 2S_{k} (-2D'_{k}) \sigma
%          
clearvars; close all;format long e;

zk = 3.0;
type = 'cgrph';

pref = [];pref.k = 16;
ns = 1;nt = 1;
ppw = 10;   % points per wavelength;
maxchunklen = pref.k/ppw/real(zk)*2*pi;

[chnkr, sources, targets] = get_geometry(type, pref, ns, nt, maxchunklen);
fprintf('Done building geometry\n');

% Plot the boundary
figure(1)
plot(chnkr, 'b.');
axis equal

% Set up kernels
Dkp = kernel('helm', 'dprime', zk);
Sk  = kernel('helm', 's', zk);
Z   = kernel.zeros();

c = 2.0;
K = [ Z       c*Sk;
      c*Dkp   Z ];
K = kernel(K);
Keval = kernel([Z c*Sk]);

% Set up boundary data
src = chnkr.r(:,:);x=src(1,:);
ubdry = 4*x.^3+2*x.^2-3*x-1;

npts = chnkr.npt
nsys = K.opdims(1)*npts;
rhs = zeros(nsys, 1);
rhs(1:K.opdims(1):end) = ubdry;

% Build the system matrix
opts = [];opts.l2scale = false;opts.rcip = true;
opts.nsub_or_tol = 30;
opts.open_arc_eye = true;
start = tic;
A = chunkermat(chnkr, K, opts) + eye(nsys);
t1 = toc(start);
fprintf('%5.2e s : time to build the system matrix\n', t1)

% Solve the linear system
start = tic;
sol = gmres(A, rhs, [], 1e-12, 200);
t1 = toc(start);

% Compute the numerical solution 
opts.forcesmooth = false;
opts.verb = false;
opts.quadkgparams = {'RelTol', 1e-8, 'AbsTol', 1.0e-8};

if isa(chnkr, 'chunkgraph')
    chnkrs = chnkr.echnks;
    chnkrtotal = merge(chnkrs);
else
    chnkrtotal = chnkr;
end

start = tic;
unum = chunkerkerneval(chnkrtotal, Keval, sol, targets, opts)
t2 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n', t2)

% Reference solution by Johan Helsing
uref = 0.02788626934981090 - 1i*0.75932847390327920
relerr  = norm(unum-uref) / norm(uref)

function [chnkobj, sources, targets] = get_geometry(type, pref, ns, nt, maxchunklen)
if strcmpi(type, 'cgrph')
    nverts = 2; 
    verts = [-1 1;-0.2 -0.2];
    edge2verts = [1;2];

    fchnks = [];
    cparams = [];
    cparams.nover = 2;
    cparams.maxchunklen = maxchunklen;
    cparams.ta = 0; cparams.tb = 1;

    chnkobj = chunkgraph(verts, edge2verts, fchnks, cparams, pref);
    chnkobj = balance(chnkobj);
    
    targets = [0.17;0.62];sources = [];
end
end

