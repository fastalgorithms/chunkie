addpaths_loc();
clear();

zk = 6*pi;

type = 'chnkr-star';
% type = 'chnkr-torus';
% type = 'cgrph';
% type = 'cgrph-sphere';

irep = 'rpcomb';

pref = [];
pref.k = 16;
ns = 10;
nt = 100;
ppw = 10;   % points per wavelength;
maxchunklen = pref.k/ppw/real(zk)*2*pi;

[chnkr, sources, targets] = get_geometry(type, pref, ns, nt, maxchunklen);
wts = chnkr.wts; wts = wts(:);

fprintf('Done building geometry\n');

% source strengths
strengths = randn(ns, 1);


% Plot everything

figure(1)
clf
hold off
plot(chnkr, 'k.');
hold on
quiver(chnkr);
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


% Set up boundary data
Skp    = kernel('axissymhelm', 'sprime', zk);
srcinfo  = []; srcinfo.r = sources; 
targinfo = []; targinfo.r = chnkr.r(:,:); targinfo.n = chnkr.n(:,:);
kernmats = Skp.eval(srcinfo, targinfo);
ubdry = kernmats*strengths;


[Skpmat, Sikmat, Sikpmat, Dkdiffmat] = get_neumann_matrices(zk, chnkr);
[Skpmat2, Sikmat2, Sikpmat2, Dkdiffmat2] = get_neumann_matrices2(zk, chnkr);


alpha = 1;
c1 = -1/(0.5 + 1i*alpha*0.25);
c2 = c1*1i*alpha;
c3 = -1;

n = chnkr.npt;
start = tic;
M = c1*(Skpmat + 1i*alpha*(Dkdiffmat*Sikmat + Sikpmat*Sikpmat)) + eye(n);
t1 = toc(start);
fprintf('Time taken for matrix matrix product:%d \n',t1);

rhs = ubdry;
start = tic; sol = M\rhs; t1 = toc(start);
fprintf('Time taken in solve: %d\n',t1);

start = tic; Minv = inv(M); t1 = toc(start);
fprintf('Time taken in computing inverse: %d\n',t1);

start = tic; sol2 = Minv*rhs; t1 = toc(start);
fprintf('Time taken in computing solution with precomp inverse: %d\n',t1);

siksol = Sikmat*sol;

%% Begin postprocessing


Sk    = kernel('axissymhelm', 's', zk);
Dk    = kernel('axissymhelm', 'd', zk);

% Compute exact solution
srcinfo  = []; srcinfo.r  = sources;
targinfo = []; targinfo.r = targets;
kernmatstarg = Sk.eval(srcinfo, targinfo);
utarg = kernmatstarg*strengths;


Keval = c1*kernel([Sk 1i*alpha*Dk]);


soluse = zeros(2*n,1);
soluse(1:2:end) = sol;
soluse(2:2:end) = siksol;

if isa(chnkr, 'chunkgraph')
    % Collapse cgrph into chnkrtotal
    chnkrs = chnkr.echnks;
    chnkrtotal = merge(chnkrs);
else
    chnkrtotal = chnkr;
end


start = tic;
Dsol = chunkerkerneval(chnkrtotal, Keval, soluse, targets);
t2 = toc(start);
fprintf('%5.2e s : time to eval at targs (slow, adaptive routine)\n', t2);

wchnkr = chnkrtotal.wts; wchnkr = wchnkr(:);
relerr  = norm(utarg-Dsol) / (sqrt(chnkrtotal.nch)*norm(utarg));
relerr2 = norm(utarg-Dsol, 'inf') / dot(abs(sol(:)), wchnkr(:));
fprintf('relative frobenius error %5.2e\n', relerr);
fprintf('relative l_inf/l_1 error %5.2e\n', relerr2);

 
 
 
function [Skpmat, Sikmat, Sikpmat, Dkdiffmat] = get_neumann_matrices(zk, chnkr)
    
    nn = chnkr.npt;
    Skp    = kernel('axissymhelm', 'sprime', zk);
    Sik    = kernel('axissymhelm', 's', 1i*zk);
    Sikp   = kernel('axissymhelm', 'sprime', 1i*zk);
    Dkdiff = kernel('axissymhelmdiff', 'dprime', [zk 1i*zk]);
    
    opts = [];
    opts.nonsmoothonly = true;
    opts.rcip = false;
    opts.nsub_or_tol = 20;
    start = tic;
    
    kernels = cell(1,4);
    kernels{1} = Skp;
    kernels{2} = Sik;
    kernels{3} = Sikp;
    kernels{4} = Dkdiff;
    
    spmats = cell(1,4);
    
    for imat=1:4
        spmats{imat} = chunkermat(chnkr, kernels{imat}, opts); 
    end
    Skp_spmat = spmats{1}; 
    Sik_spmat = spmats{2}; 
    Sikp_spmat = spmats{3};
    Dkdiff_spmat = spmats{4};


    
    t1 = toc(start);
    fprintf('Time taken in sparse matrix generation: %d\n',t1);
    
    nthd = maxNumCompThreads;
    nbsize = ceil(chnkr.npt/nthd);
    Skp_cellmat = cell(nthd,1);
    Sik_cellmat = cell(nthd,1);
    Sikp_cellmat =  cell(nthd,1);
    Dkdiff_cellmat = cell(nthd,1);
    l2scale = false;
    start = tic;
    nn = chnkr.npt;
    for i=1:nthd 
        opdims = [1 1];
        iind = ((i-1)*nbsize+1):min(nn,i*nbsize);
        Skp_cellmat{i} = chnk.flam.kernbyindex(iind, 1:nn, chnkr, ... 
                Skp, opdims, Skp_spmat);
            
            
        Sik_cellmat{i} = chnk.flam.kernbyindex(iind, 1:nn, chnkr, ... 
                Sik, opdims, Sik_spmat);
            
            
        Sikp_cellmat{i} = chnk.flam.kernbyindex(iind, 1:nn, chnkr, ... 
                Sikp, opdims, Sikp_spmat);
            
            
        Dkdiff_cellmat{i} = chnk.flam.kernbyindex(iind, 1:nn, chnkr, ... 
                Dkdiff, opdims, Dkdiff_spmat);
    end
    
    Skpmat = vertcat(Skp_cellmat{:});
    Sikmat = vertcat(Sik_cellmat{:});
    Sikpmat = vertcat(Sikp_cellmat{:});
    Dkdiffmat = vertcat(Dkdiff_cellmat{:});
    t1 = toc(start);
    fprintf('Time taken in matrix entry generation: %d\n',t1);
end


function [Skpmat, Sikmat, Sikpmat, Dkdiffmat] = get_neumann_matrices2(zk, chnkr)
    
 
    Nkerns = kernel('axissymhelm', 'neu_rpcomb', zk);
    
    opts = [];
    opts.nonsmoothonly = true;
    opts.rcip = false;
    opts.nsub_or_tol = 20;
    start = tic;
    
    spmat = chunkermat(chnkr, Nkerns, opts);
    
    
    t1 = toc(start);
    fprintf('Time taken in sparse matrix generation (new): %d\n',t1);
    
    nsys = 3*chnkr.npt;
    nthd = maxNumCompThreads;
    nbsize = ceil(nsys/nthd);
    N_cellmat = cell(nthd,1);
    start = tic;
    for i=1:nthd 
        opdims = Nkerns.opdims;
        iind = ((i-1)*nbsize+1):min(nsys,i*nbsize);
        N_cellmat{i} = chnk.flam.kernbyindex(iind, 1:nsys, chnkr, ... 
                Nkerns, opdims, spmat);           
    end
    
    Nmat = vertcat(N_cellmat{:});
    Skpmat = Nmat(1:3:end, 1:3:end)/Nkerns.params.c1;
    Dkdiffmat = Nmat(1:3:end, 2:3:end)/Nkerns.params.c2;
    Sikpmat = Nmat(1:3:end, 3:3:end)/Nkerns.params.c2;
    Sikmat = -Nmat(2:3:end, 1:3:end);
    t1 = toc(start);
    fprintf('Time taken in matrix entry generation (new): %d\n',t1);
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
    cparams.ta = -pi/2;
    cparams.tb = pi/2;
    cparams.maxchunklen = maxchunklen;
    narms = 0;
    amp = 0.0;
    
    fchnks{1} = @(t) circle(t);
    
    chnkobj = chunkgraph(verts, edge2verts, fchnks, cparams, pref);
    chnkobj = balance(chnkobj);
    
    ts = -pi/2 + pi*rand(ns, 1);
    sources = starfish(ts, narms, amp);
    sources = 0.1*sources;

    ts = -pi/2 + pi*rand(nt, 1);
    targets = starfish(ts, narms, amp);
    targets = targets .* (1 + 0.5*repmat(rand(1, nt), 2, 1));   

    
    
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
xs = cos(-pi/4 + pi/2*t);
ys = sin(-pi/4 + pi/2*t);
xp = -pi*ys/2;
yp = pi*xs/2;
xpp = -pi*pi*xs/4;
ypp = -pi*pi*ys/4;
r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end



