pquadTest0();
pquadTest1();
pquadTest2();
pquadTest3();
pquadTest4();

function pquadTest0()
% testing product quadrature rule
% planewave vec

kvec = 10*[1;-1.5];

zk = norm(kvec);

% define geometry and boundary conditions
% (vertices defined by another function)

cparams = []; cparams.eps = 1e-9;
chnkr = chunkerfunc(@(t) starfish(t),cparams);
chnkr = refine(chnkr,struct('nover',1));

% solve and visualize the solution

% build laplace dirichlet matrix

fkern = kernel('helmholtz','c',zk,[1.5,-2i]);
opts = [];
start = tic; C = chunkermat(chnkr,fkern,opts);
t1 = toc(start);

fprintf('%5.2e s : time to assemble matrix\n',t1)

opts.forcepquad = true;
opts.side = 'i';
sys = chunkerkernevalmat(chnkr,fkern,[chnkr.r(1,:);chnkr.r(2,:)],opts);

rhs = besselh(0,zk*abs((1+1.25*1i)-(chnkr.r(1,:)+1i*chnkr.r(2,:)))); rhs = rhs(:);
start = tic; sol = gmres(sys,rhs,[],1e-14,100); t1 = toc(start);

fprintf('%5.2e s : time for dense gmres\n',t1)

% evaluate at targets and plot

rmin = min(chnkr); rmax = max(chnkr);
nplot = 300;
hx = (rmax(1)-rmin(1))/nplot;
hy = (rmax(2)-rmin(2))/nplot;
xtarg = linspace(rmin(1)+hx/2,rmax(1)-hx/2,nplot);
ytarg = linspace(rmin(2)+hy/2,rmax(2)-hy/2,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);

start = tic; in = chunkerinterior(chnkr,targets); t1 = toc(start);
fprintf('%5.2e s : time to find points in domain\n',t1)

% compute layer potential at interior points
opts.forcepquad = true;
opts.side = 'i'; % 'i' for interior, 'e' for exterior, for positively oriented curve.
start = tic;
Csolpquad = chunkerkerneval(chnkr,fkern,sol,targets(:,in),opts);
t1 = toc(start);
fprintf('%5.2e s : time for kerneval (Helsing-Ojala for near)\n',t1);

start = tic;
Csol = chunkerkerneval(chnkr,fkern,sol,targets(:,in)); t1 = toc(start);
fprintf('%5.2e s : time for kerneval (adaptive for near)\n',t1);

% Compare with reference solution Dsol
rel_error = max(abs(Csol-Csolpquad))/max(abs(Csol));
fprintf('%5.2e : Relative max error\n',rel_error);

assert(rel_error < 1e-10)

% same combined kernel, built by arithmetic, should agree under pquad
fkernd = kernel('helmholtz','d',zk);
fkerns = kernel('helmholtz','s',zk);
fkernc_arith = 1.5*fkernd + (-2i)*fkerns;
assert(~isempty(fkernc_arith.splitinfo))

Csolpquad_arith = chunkerkerneval(chnkr,fkernc_arith,sol,targets(:,in),opts);
rel_error_arith = max(abs(Csolpquad-Csolpquad_arith))/max(abs(Csolpquad));
fprintf('%5.2e : Relative max error, arithmetic vs native combined kernel\n',rel_error_arith);

assert(rel_error_arith < 1e-10)

end

function pquadTest1()
% check that scalar times and minus also combine splitinfo correctly

zk = 1.3;
coefs = [1.5,-2i];

cparams = []; cparams.eps = 1e-9;
chnkr = chunkerfunc(@(t) starfish(t),cparams);

fkernc = kernel('helmholtz','c',zk,coefs);
fkernd = kernel('helmholtz','d',zk);
fkerns = kernel('helmholtz','s',zk);

fkernc_arith = coefs(1)*fkernd - (-coefs(2))*fkerns;
assert(~isempty(fkernc_arith.splitinfo))

opts = [];
opts.forcepquad = true;
opts.side = 'i';
targs = chnkr.r(:,:) - 1e-3*chnkr.n(:,:);

sysc = chunkerkernevalmat(chnkr,fkernc,targs,opts);
sysc_arith = chunkerkernevalmat(chnkr,fkernc_arith,targs,opts);

rel_error = norm(sysc-sysc_arith,'fro')/norm(sysc,'fro');
fprintf('%5.2e : Relative fro error in native vs arithmetic\n',rel_error);

assert(rel_error < 1e-10)

end

function pquadTest2()
% check that kernel([a,b]) and kernel([a;b]) reproduce block row/column
% kernel evaluation

zk = 0.9;
a = kernel('lap','s');
b = kernel('helm','d',zk);

src = []; src.r = [[0;0],[0.4;1]];
src.d = randn(2,2); src.d2 = randn(2,2); src.n = randn(2,2);
targ = []; targ.r = [[2;2],[-1;0.4],[0.1;-0.7]];
targ.d = randn(2,3); targ.d2 = randn(2,3); targ.n = randn(2,3);

Arow = a.eval(src,targ); Brow = b.eval(src,targ);
ref_row = zeros(size(targ.r,2),2*size(src.r,2));
ref_row(:,1:2:end) = Arow; ref_row(:,2:2:end) = Brow;

Krow = kernel([a,b]);
assert(norm(Krow.eval(src,targ)-ref_row) < 1e-10)

ref_col = zeros(2*size(targ.r,2),size(src.r,2));
ref_col(1:2:end,:) = Arow; ref_col(2:2:end,:) = Brow;

Kcol = kernel([a;b]);
assert(norm(Kcol.eval(src,targ)-ref_col) < 1e-10)

end

function pquadTest3()
% check mtimes splitinfo: pointwise scalar (opdims=[1,1]), and left/right
% by non-square matrices on a matrix-valued (opdims=[2,2]) kernel

zk = 1.1;
kernd = kernel('helmholtz','d',zk);
mu = 1.2;
kernv = kernel('stokes','svel',mu);

cparams = []; cparams.eps = 1e-9;
chnkr = chunkerfunc(@(t) starfish(t),cparams);
targs = chnkr.r(:,:) - 1e-3*chnkr.n(:,:);
nt = size(targs,2);
ns = chnkr.npt;

opts = [];
opts.forcepquad = true;
opts.side = 'i';

% pointwise scalar function, left multiply: splitinfo is preserved
wfun = @(t) reshape(1 + t.r(1,:).^2, 1, 1, []);
wkernd = wfun*kernd;
assert(~isempty(wkernd.splitinfo))

sysd = chunkerkernevalmat(chnkr,kernd,targs,opts);
syswd = chunkerkernevalmat(chnkr,wkernd,targs,opts);

w = (1 + targs(1,:).^2).';
rel_error_w = norm(w.*sysd - syswd,'fro')/norm(syswd,'fro');
fprintf('%5.2e : Relative fro error in w(t)*K\n',rel_error_w);
assert(rel_error_w < 1e-10)

% right multiply by a function handle: splitinfo is dropped, since
% pquad can't handle a split coefficient that varies within a panel
wrkernd = kernd*wfun;
assert(isempty(wrkernd.splitinfo))

% matrix-valued kernel, left multiply by a 3x2 constant matrix: p~=m
sysv = chunkerkernevalmat(chnkr,kernv,targs,opts);

A = [1, 2; 3, -1; 0, 1];
Akernv = A*kernv;
assert(~isempty(Akernv.splitinfo))
sysA = chunkerkernevalmat(chnkr,Akernv,targs,opts);

rel_error_left = norm(kron(eye(nt),A)*sysv - sysA,'fro')/norm(sysA,'fro');
fprintf('%5.2e : Relative fro error in A*K\n',rel_error_left);
assert(rel_error_left < 1e-10)

% matrix-valued kernel, right multiply by a 2x3 constant matrix: p~=q
B = [1, 0, -2; 2, 1, 1];
kernvB = kernv*B;
assert(~isempty(kernvB.splitinfo))
sysB = chunkerkernevalmat(chnkr,kernvB,targs,opts);

rel_error_right = norm(sysv*kron(eye(ns),B) - sysB,'fro')/norm(sysB,'fro');
fprintf('%5.2e : Relative fro error in K*B\n',rel_error_right);
assert(rel_error_right < 1e-10)

end

function pquadTest4()
% check that the corrections=true option is correct

zk = 1.7;
fkern = kernel('helm','s',zk);

cparams = []; cparams.eps = 1e-9;
chnkr = chunkerfunc(@(t) starfish(t),cparams);
chnkr = refine(chnkr,struct('nover',1));

targinfo = []; targinfo.r = chnkr.r(1:2,:) - 1e-3*chnkr.n(1:2,:);

opts = [];
opts.forcepquad = true;
opts.side = 'i';

full_fp = chunkerkernevalmat(chnkr,fkern,targinfo.r,opts);

opts_c = opts; opts_c.corrections = true;
cor_fp = chunkerkernevalmat(chnkr,fkern,targinfo.r,opts_c);

srcinfo = []; srcinfo.r = chnkr.r(:,:); srcinfo.n = chnkr.n(:,:);
srcinfo.d = chnkr.d(:,:); srcinfo.d2 = chnkr.d2(:,:);
wts = chnkr.wts; wts = wts(:).';
smooth_all = fkern.eval(srcinfo,targinfo).*wts;

rel_error = norm(full_fp - (smooth_all+cor_fp),'fro')/norm(full_fp,'fro');
fprintf('%5.2e : Relative fro error, full vs smooth + corrections (forcepquad)\n',rel_error);

assert(rel_error < 1e-10)

end
