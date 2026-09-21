kernel_onesTest0();
kernel_onesTest1();
kernel_onesTest2();


function kernel_onesTest0()
%KERNEL_ONESTEST0 test kernel.ones: constructor, eval, and fmm.

iseed = 8675309;
rng(iseed);

% basic properties

K = kernel('ones');
assert(isequal(K.opdims, [1 1]));
A = [1 2 3; 4 5 6];
KA = kernel.ones(A);
assert(isequal(KA.opdims, [2 3]));
s = []; s.r = randn(2, 4); t = []; t.r = randn(2, 5);
assert(isequal(KA.eval(s, t), repmat(A, 5, 4)));
sig = randn(3*4, 1);
assert(norm(KA.fmm(1e-10, s, t, sig) - KA.eval(s, t)*sig) < 1e-12);

% fmm with a single output and multiple outputs (grad/hess are zero)
K1 = kernel.ones();
sig = randn(4, 1);
pot = K1.fmm(1e-10, s, t, sig);
assert(norm(pot - sum(sig)*ones(5, 1)) < 1e-12);
[pot, grad, hess] = K1.fmm(1e-10, s, t, sig);
assert(isequal(size(grad), [2 5]) && all(grad(:) == 0));
assert(isequal(size(hess), [3 5]) && all(hess(:) == 0));

% fmm with targets given as a plain array
pot2 = KA.fmm(1e-10, s, t.r, randn(3*4, 1));
assert(isequal(size(pot2), [2*5 1]));

% shifted_eval agrees with eval
assert(isequal(KA.shifted_eval(s, t, [0.3; -0.2]), KA.eval(s, t)));

end


function kernel_onesTest1()
%KERNEL_ONESTEST1 interior Laplace Neumann solve where the ones kernel
% removes the null space of 1/2 I + S'.

iseed = 8675309;
rng(iseed);

% geometry

cparams = []; cparams.eps = 1.0e-10;
pref = []; pref.k = 16;
narms = 3; amp = 0.25;
chnkr = chunkerfunc(@(t) starfish(t, narms, amp), cparams, pref);

% exact harmonic solution from a unit point charge outside the curve

src = [2.5; 1.3];

kerns = kernel('lap', 's');
srcinfo = []; srcinfo.r = src;
kernsp = kernel('lap', 'sprime');
rhs = kernsp.eval(srcinfo, chnkr);

% Add the rank-one ones kernel W sigma = int sigma to remove the null space.

sysmat = chunkermat(chnkr, kernel('lap', 'sprime') + kernel('ones'));
sysmat = sysmat + 0.5*eye(chnkr.npt);

assert(cond(sysmat) < 1e3, 'ones kernel failed to remove the null space');

sol = sysmat \ rhs;

% the ones term forces int sigma = 0
assert(abs(sum(chnkr.wts(:).*sol)) < 1e-8);

% random targets in the box 
nt = 20;
targs = []; targs.r = rand(2, nt) - 0.5;

uex = kerns.eval(srcinfo, targs);
u = chunkerkerneval(chnkr, kernel('lap', 's'), sol, targs);

% Neumann solution determined up to a constant
uex = uex - mean(uex); u = u - mean(u);
relerr = norm(u - uex)/norm(uex);
fprintf('relative error interior Neumann (ones kernel): %5.2e\n', relerr);
assert(relerr < 1e-8);

% ones kernel with fmm path in chunkerkerneval
w = chunkerkerneval(chnkr, kernel('ones'), sol, targs);
assert(norm(w) < 1e-8);

end


function kernel_onesTest2()
%KERNEL_ONESTEST2 exterior Laplace Dirichlet solve, using the representation
%
%   u = D[sigma] + G(x,0) * int sigma
%
% The ones kernel, scaled by G(t,0), supplies the log term at infinity

iseed = 8675309;
rng(iseed);

% geometry

chnkr = chunkerfunc(@(t) starfish(t));

src0 = []; src0.r = [0; 0];

% exact solution from a unit point charge inside the curve (nonzero
% total charge, so the log term is required)

srcs = []; srcs.r = [0.3; -0.1];
strengths = 1.0;

skern = kernel('lap', 's');
rhs = skern.eval(srcs, chnkr)*strengths;

Gt = @(t) reshape(skern.eval(src0, t), 1, 1, []);
kerns = kernel('lap', 'd') + Gt * kernel.ones();

sysmat = chunkermat(chnkr, kerns);
sysmat = sysmat + 0.5*eye(chnkr.npt);

sigma = sysmat \ rhs;

% evaluate at exterior targets

targs = [[1.7; 2.1], [-2.5; 0.4], [0.3; -3.0]];
tinfo = []; tinfo.r = targs;
pot_ex = skern.eval(srcs, tinfo)*strengths;

opts = []; opts.forcefmm = true;
pot = chunkerkerneval(chnkr, kerns, sigma, targs, opts);
err = norm(pot - pot_ex)/norm(pot_ex);
fprintf('Error in exterior Dirichlet solve = %5.2e\n', err);
assert(err < 1e-8, 'exterior Dirichlet error too large');

end

