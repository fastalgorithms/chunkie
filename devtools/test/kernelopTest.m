%KERNELOPTEST verify kernel class operations are correct
kernelopTest0()
kernelopTest1()
kernelopTest2()
kernelopTest3()

function kernelopTest0()
% setup points
src = [];
src.r = [[0;0],[0;1]];
src.n = randn(size(src.r)); src.n = src.n./vecnorm(src.n);

targ = [];
targ.r = [[1;1],[-1;0]];
targ.n = randn(size(targ.r)); targ.n = targ.n./vecnorm(targ.n);

% set primitive kernels
skern = kernel('lap','s');
dkern = kernel('helm','d',1);

% set scalar
a = pi*1i;

% interleave kernels
fkern1 =kernel([skern;dkern]);

% multiply and divide by a
fkern2 = a.*fkern1;
fkern3 = fkern1/a;

assert(norm(fkern2.eval(src,targ) - a*fkern1.eval(src,targ))<1e-10)
assert(norm(fkern3.eval(src,targ) - 1/a*fkern1.eval(src,targ))<1e-10)

% negate a kernel
nkern = -skern;
assert(norm(nkern.eval(src,targ) + skern.eval(src,targ))<1e-10)

% add and subtract kernels
ckern1 = skern + dkern;
ckern2 = skern - dkern;
assert(norm(ckern1.eval(src,targ) - (skern.eval(src,targ)+dkern.eval(src,targ)))<1e-10)
assert(norm(ckern2.eval(src,targ) - (skern.eval(src,targ)-dkern.eval(src,targ)))<1e-10)

% conjugate a kernel
conj_dkern = conj(dkern);
assert(norm(conj_dkern.eval(src,targ) - conj_dkern.eval(src,targ))<1e-10)

end

function kernelopTest1()
% Test matrix-valued mtimes for a vector kernel (Stokes velocity,
% opdims = [2 2]): left and right multiplication by scalar and
% matrix-valued functions and constants, square and non-square,
% checked against eval and fmm.

% setup points
src = [];
src.r = [[0;0],[0;1]];
src.n = randn(size(src.r)); src.n = src.n./vecnorm(src.n);

targ = [];
targ.r = [[1;1],[-1;0]];
targ.n = randn(size(targ.r)); targ.n = targ.n./vecnorm(targ.n);

stkern = kernel('stokes','svel',1.0);   % opdims = [2 2]
Kmat = stkern.eval(src, targ);   % (2*nt x 2*ns)
nt = size(targ.r, 2);
ns = size(src.r, 2);

% larger point set, for fmm checks
rng(1);
nsrc = 12; ntarg = 9;
src2 = []; src2.r = randn(2, nsrc);
src2.n = randn(2, nsrc); src2.n = src2.n./vecnorm(src2.n);
targ2 = []; targ2.r = 3 + randn(2, ntarg);
targ2.n = randn(2, ntarg); targ2.n = targ2.n./vecnorm(targ2.n);

sigma = randn(2*nsrc, 1);
eps_fmm = 1e-12;
Kfmm_out = stkern.fmm(eps_fmm, src2, targ2, sigma);   % (2*ntarg x 1)

%% scalar constant: left and right

scaled_l = 2 * stkern;
assert(norm(scaled_l.eval(src,targ) - 2*Kmat) < 1e-10, ...
    'scalar*K mtimes mismatch');

scaled_r = stkern * 2;
assert(norm(scaled_r.eval(src,targ) - 2*Kmat) < 1e-10, ...
    'K*scalar mtimes mismatch');

%% scalar function (pointwise multiplier): left and right

% w(t) = 1 + t.r(1,:).^2, a scalar function of the target
wfun_t = @(t) reshape(1 + t.r(1,:).^2, 1, 1, []);
wKfun = wfun_t * stkern;
assert(isequal(wKfun.opdims, stkern.opdims), 'w(t)*K should preserve opdims');

% eval
ref_w_left = zeros(2*nt, 2*ns);
for j = 1:nt
    wj = 1 + targ.r(1,j)^2;
    ref_w_left((j-1)*2+(1:2), :) = wj * Kmat((j-1)*2+(1:2), :);
end
got_w_left = wKfun.eval(src, targ);
assert(norm(got_w_left - ref_w_left) < 1e-10, 'w(t)*K eval mismatch');

% fmm
wKfun_fmm = wKfun.fmm(eps_fmm, src2, targ2, sigma);
ref_wfmm_left = zeros(2*ntarg, 1);
for j = 1:ntarg
    wj = 1 + targ2.r(1,j)^2;
    ref_wfmm_left((j-1)*2+(1:2)) = wj * Kfmm_out((j-1)*2+(1:2));
end
assert(norm(wKfun_fmm - ref_wfmm_left) < 1e-8, 'pointwise fmm_left mismatch');

% w(s) = 1 + s.r(2,:).^2, a scalar function of the source
wfun_s = @(s) reshape(1 + s.r(2,:).^2, 1, 1, []);
Kwfun = stkern * wfun_s;
assert(isequal(Kwfun.opdims, stkern.opdims), 'K*w(s) should preserve opdims');

% eval
ref_w_right = zeros(2*nt, 2*ns);
for i = 1:ns
    wi = 1 + src.r(2,i)^2;
    ref_w_right(:, (i-1)*2+(1:2)) = Kmat(:, (i-1)*2+(1:2)) * wi;
end
got_w_right = Kwfun.eval(src, targ);
assert(norm(got_w_right - ref_w_right) < 1e-10, 'K*w(s) eval mismatch');

% fmm
Kwfun_fmm = Kwfun.fmm(eps_fmm, src2, targ2, sigma);
sigma3 = zeros(2*nsrc, 1);
for i = 1:nsrc
    wi = 1 + src2.r(2,i)^2;
    sigma3((i-1)*2+(1:2)) = wi * sigma((i-1)*2+(1:2));
end
ref_wfmm_right = stkern.fmm(eps_fmm, src2, targ2, sigma3);
assert(norm(Kwfun_fmm - ref_wfmm_right) < 1e-8, 'pointwise fmm_right mismatch');

%% square matrix constant: left and right

A = [1, 2; 3, 4];

AKern = A * stkern;
assert(isequal(AKern.opdims, [2, 2]));
assert(norm(AKern.eval(src, targ) - kron(eye(nt), A)*Kmat) < 1e-10, ...
    'A*K eval mismatch');

KAkern = stkern * A;
assert(isequal(KAkern.opdims, [2, 2]));
assert(norm(KAkern.eval(src, targ) - Kmat*kron(eye(ns), A)) < 1e-10, ...
    'K*A eval mismatch');

%% square matrix-valued function: left and right

% Left multiply by a 2x2 matrix-valued function M(t) = [1, 0; -t.r(1,:), 1].
Mfun = @(t) reshape([ones(1,size(t.r,2)); -t.r(1,:); ...
                     zeros(1,size(t.r,2)); ones(1,size(t.r,2))], 2, 2, []);
MKfun = Mfun * stkern;
assert(isequal(MKfun.opdims, [2, 2]));

% eval
ref_left = zeros(2*nt, 2*ns);
for j = 1:nt
    kappa_j = targ.r(1,j);
    Mj = [1, 0; -kappa_j, 1];
    ref_left((j-1)*2+(1:2), :) = Mj * Kmat((j-1)*2+(1:2), :);
end
got_left = MKfun.eval(src, targ);
assert(norm(got_left - ref_left) < 1e-10, 'M(t)*K eval mismatch');

% fmm
MKfun_fmm = MKfun.fmm(eps_fmm, src2, targ2, sigma);
ref_fmm_left = zeros(2*ntarg, 1);
for j = 1:ntarg
    kappa_j = targ2.r(1,j);
    Mj = [1, 0; -kappa_j, 1];
    ref_fmm_left((j-1)*2+(1:2)) = Mj * Kfmm_out((j-1)*2+(1:2));
end
assert(norm(MKfun_fmm - ref_fmm_left) < 1e-8, 'fmm_left mismatch');

% Right multiply by a 2x2 matrix-valued function N(s) = [1, 0; -s.r(2), 1].
Nfun = @(s) reshape([ones(1,size(s.r,2)); -s.r(2,:); ...
                     zeros(1,size(s.r,2)); ones(1,size(s.r,2))], 2, 2, []);
KNfun = stkern * Nfun;
assert(isequal(KNfun.opdims, [2, 2]));

% eval
ref_right = zeros(2*nt, 2*ns);
for i = 1:ns
    kappa_i = src.r(2,i);
    Ni = [1, 0; -kappa_i, 1];
    ref_right(:, (i-1)*2+(1:2)) = Kmat(:, (i-1)*2+(1:2)) * Ni;
end
got_right = KNfun.eval(src, targ);
assert(norm(got_right - ref_right) < 1e-10, 'K*N(s) eval mismatch');

% fmm
KNfun_fmm = KNfun.fmm(eps_fmm, src2, targ2, sigma);
sigma2 = zeros(2*nsrc, 1);
for i = 1:nsrc
    kappa_i = src2.r(2,i);
    Ni = [1, 0; -kappa_i, 1];
    sigma2((i-1)*2+(1:2)) = Ni * sigma((i-1)*2+(1:2));
end
ref_fmm_right = stkern.fmm(eps_fmm, src2, targ2, sigma2);
assert(norm(KNfun_fmm - ref_fmm_right) < 1e-8, 'fmm_right mismatch');

%% non-square matrix constant: left and right

% B is 1x2
B = [1, -2];
BKern = B * stkern;
assert(isequal(BKern.opdims, [1, 2]));
assert(norm(BKern.eval(src, targ) - kron(eye(nt), B)*Kmat) < 1e-10, ...
    'B*K (non-square constant) eval mismatch');

% C is 2x1
C = [1; -2];
KCkern = stkern * C;
assert(isequal(KCkern.opdims, [2, 1]));
assert(norm(KCkern.eval(src, targ) - Kmat*kron(eye(ns), C)) < 1e-10, ...
    'K*C (non-square constant) eval mismatch');

%% non-square matrix-valued function: left and right

% P(t) is 3x2
Pfun = @(t) reshape([ones(1,size(t.r,2)); -t.r(1,:); t.r(2,:); ...
                     zeros(1,size(t.r,2)); ones(1,size(t.r,2)); ...
                     2*ones(1,size(t.r,2))], 3, 2, []);
PKfun = Pfun * stkern;
assert(isequal(PKfun.opdims, [3, 2]));

ref_P_left = zeros(3*nt, 2*ns);
for j = 1:nt
    Pj = [1, 0; -targ.r(1,j), 1; targ.r(2,j), 2];
    ref_P_left((j-1)*3+(1:3), :) = Pj * Kmat((j-1)*2+(1:2), :);
end
got_P_left = PKfun.eval(src, targ);
assert(norm(got_P_left - ref_P_left) < 1e-10, 'P(t)*K (non-square) eval mismatch');

% Q(s) is 2x3
Qfun = @(s) reshape([ones(1,size(s.r,2)); -s.r(2,:); ...
                     zeros(1,size(s.r,2)); ones(1,size(s.r,2)); ...
                     s.r(1,:); 2*ones(1,size(s.r,2))], 2, 3, []);
KQfun = stkern * Qfun;
assert(isequal(KQfun.opdims, [2, 3]));

ref_Q_right = zeros(2*nt, 3*ns);
for i = 1:ns
    Qi = [1, 0, src.r(1,i); -src.r(2,i), 1, 2];
    ref_Q_right(:, (i-1)*3+(1:3)) = Kmat(:, (i-1)*2+(1:2)) * Qi;
end
got_Q_right = KQfun.eval(src, targ);
assert(norm(got_Q_right - ref_Q_right) < 1e-10, 'K*Q(s) (non-square) eval mismatch');

end

function kernelopTest2()
% Test mtimes shifted_eval (axisymmetric Helmholtz kernel, opdims = [1 1]):
% left and right multiplication by scalar functions.

zk = 1.3;
axkern = kernel('axissymhelm', 's', zk);   % opdims = [1 1], has shifted_eval

% axisymmetric kernel requires nonnegative radial (first) coordinate
asrc = []; asrc.r = [[0.5;0],[1.0;1]];
atarg = []; atarg.r = [[2.0;1],[1.5;-0.5]];
nat = size(atarg.r, 2); nas = size(asrc.r, 2);

o = [0.7, -0.3];
Kshift = axkern.shifted_eval(asrc, atarg, o);

% pointwise scalar multiplier, left
wfun_t2 = @(t) reshape(1 + t.r(1,:).^2, 1, 1, []);
axMK = wfun_t2 * axkern;
assert(isequal(axMK.opdims, axkern.opdims));

ref_shift_left = zeros(nat, nas);
for j = 1:nat
    wj = 1 + (atarg.r(1,j) + o(1))^2;
    ref_shift_left(j, :) = wj * Kshift(j, :);
end
got_shift_left = axMK.shifted_eval(asrc, atarg, o);
assert(norm(got_shift_left - ref_shift_left) < 1e-10, 'shifted_eval_left mismatch');

% pointwise scalar multiplier, right
wfun_s2 = @(s) reshape(1 + s.r(2,:).^2, 1, 1, []);
axKN = axkern * wfun_s2;
assert(isequal(axKN.opdims, axkern.opdims));

ref_shift_right = zeros(nat, nas);
for i = 1:nas
    wi = 1 + (asrc.r(2,i) + o(2))^2;
    ref_shift_right(:, i) = Kshift(:, i) * wi;
end
got_shift_right = axKN.shifted_eval(asrc, atarg, o);
assert(norm(got_shift_right - ref_shift_right) < 1e-10, 'shifted_eval_right mismatch');

end

function kernelopTest3()
% Test that kernel.mtimes errors on bad inputs: function handles of more
% than one argument, and matrix-valued constants/functions whose
% dimensions are inconsistent with K.opdims.

stkern = kernel('stokes','svel',1.0);   % opdims = [2 2]

% Function handles of more than one argument should be rejected
badfun = @(s,t) ones(2,2,size(s.r,2));
try
    badfun * stkern;
    error('expected error for two-argument function handle not thrown');
catch ME
    assert(strcmp(ME.identifier, 'MATLAB:assertion:failed') || ...
           contains(ME.message, 'function of source or target'), ...
        'wrong error for two-argument function handle');
end

% Constant matrix with wrong number of columns for left multiply
try
    badA = [1, 2, 3];
    badA * stkern;
    error('expected error for bad constant matrix opdims not thrown');
catch ME
    assert(contains(ME.message, 'must have'), ...
        'wrong error for bad constant matrix opdims');
end

% Function handle returning matrices with wrong inner dimension
try
    badAfun = @(t) ones(2, 3, size(t.r,2));
    badAfun * stkern;
    error('expected error for bad function handle opdims not thrown');
catch ME
    assert(contains(ME.message, 'must return matrices with'), ...
        'wrong error for bad function handle opdims');
end

end
