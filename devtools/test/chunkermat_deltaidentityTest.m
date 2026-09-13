chunkermat_deltaidentityTest0();
chunkermat_deltaidentityTest1();
chunkermat_deltaidentityTest2();


function chunkermat_deltaidentityTest0()
% kernel.eye on a smooth curve: dense, block, and kernel arithmetic

K = kernel('lap','d');
chnkr = chunkerfunc(@(t) starfish(t,3,0.25));
n = chnkr.npt;
A = chunkermat(chnkr, K);

% a*I + A
a = -0.5;
Asm = chunkermat(chnkr, K + kernel.eye(a));
err = norm(a*eye(n) + A - Asm, 'fro')/norm(A, 'fro');
assert(err < 1e-12, 'dense mismatch (err=%.2e)', err);

% 2x2 block
Z = kernel.zeros(1,1);
a = 0.5; c = 0.3;
Kd = kernel([K + kernel.eye(a), Z;...
    kernel.eye(c), 2*(K + kernel.eye(a))]);
Ad = chunkermat(chnkr, Kd);
Aref = kron(A, diag([1,2])) + kron(speye(n), [a 0; c 2*a]);
err = norm(Ad - Aref, 'fro')/norm(Ad, 'fro');
assert(err < 1e-12, 'block mismatch (err=%.2e)', err);

% arithmetic carries the delta term
a = -0.37;
Kd = K + kernel.eye(a);
D0 = Kd.diag(chnkr);

Kneg = -Kd;
assert(isequal(Kneg.diag(chnkr), -D0), 'uminus dropped diag');
assert(isequal(conj(Kd).diag(chnkr), conj(D0)), 'conj dropped diag');
Kscl = 3*Kd;
assert(isequal(Kscl.diag(chnkr), 3*D0), 'times dropped diag');
Kdiv = Kd/4;
assert(isequal(Kdiv.diag(chnkr), D0/4), 'rdivide dropped diag');
Ksub = Kd - kernel.eye(a);
assert(nnz(Ksub.diag(chnkr)) == 0, 'minus dropped diag');

% block multiply transforms the delta
M = [0 1; 2 0];
K2 = kernel([Kd, kernel.zeros(1,1); kernel.zeros(1,1), Kd]);
K2m = M*K2;
Dm = K2m.diag(chnkr);
D2 = K2.diag(chnkr);
Dref = zeros(size(D2));
for ii = 1:chnkr.npt
    Dref(2*ii-1:2*ii,:) = M*D2(2*ii-1:2*ii,:);
end
assert(isequal(Dm, Dref), 'mtimes did not transform diag');

end


function chunkermat_deltaidentityTest1()
% kernel.eye on a Laplace triangle with corners

skern = kernel('lap','s');
K = kernel('lap','d');
src.r = [3.1; 1.7];
targ.r = [0.05; -0.1];
uex = skern.eval(src, targ);

verts = [cos(2*pi*(0:2)/3 + pi/2); sin(2*pi*(0:2)/3 + pi/2)];
cg = chunkgraph(verts, [1 2 3; 2 3 1]);
rhs = skern.eval(src, cg);

% delta path matches the legacy path at a = 1
Aleg = chunkermat(cg, K);
Aone = chunkermat(cg, K + kernel.eye(1));
err = norm(Aleg + eye(size(Aleg,1)) - Aone, 'fro')/norm(Aone, 'fro');
assert(err < 1e-12, 'legacy mismatch (err=%.2e)', err);

% interior Dirichlet with a = -1/2
Khalf = K + kernel.eye(-0.5);
Ahalf = chunkermat(cg, Khalf);
sigma = Ahalf \ rhs;
u = chunkerkerneval(cg, K, sigma, targ);
err = abs(u - uex)/abs(uex);
assert(err < 1e-9, 'Dirichlet field error %.2e', err);

% complex scaling
c = -3.5i;
Ascl = chunkermat(cg, c*Khalf);
err = norm(Ascl - c*Ahalf, 'fro')/(abs(c)*norm(Ahalf, 'fro'));
assert(err < 1e-13, 'scaling mismatch (err=%.2e)', err);

% kernel.eye(0) leaves the matrix unchanged, first kind solve
A1 = chunkermat(cg, skern);
A2 = chunkermat(cg, skern + kernel.eye(zeros(1,1)));
err = norm(A1 - A2, 'fro')/norm(A1, 'fro');
assert(err < 1e-14, 'eye(0) changed the matrix (err=%.2e)', err);

sigma = A2 \ rhs;
u = chunkerkerneval(cg, skern, sigma, targ);
err = abs(u - uex)/abs(uex);
assert(err < 1e-6, 'first kind field error %.2e', err);

end


function chunkermat_deltaidentityTest2()
% Helmholtz first kind and mixed BVPs with RCIP, u = S[sigma]

zk = [2.0, 1.3];
l0 = 0.25;

nv = 3;
R = 8*l0/(2*sin(pi/nv));
ang = 2*pi*(0:nv-1)/nv;
cgh = chunkgraph(R*[cos(ang); sin(ang)], [1:nv; [2:nv, 1]], []);
cgh = refine(cgh, struct('last_len', l0));

opts = []; opts.rcip = true;

srcinfo = struct('r', [4.2; 1.4]);
targinfo = struct('r', [0.05; -0.1]);
targh = targinfo.r;

Sh = kernel('helm', 's', zk(1));
ub = chnk.helm2d.kern(zk(1), srcinfo, cgh, 's');
uexh = chnk.helm2d.kern(zk(1), srcinfo, targinfo, 's');

% Dirichlet, first kind
Ah = chunkermat(cgh, Sh + kernel.eye(0), opts);
sig = Ah \ ub;
u_t = chunkerkerneval(cgh, Sh, sig, targh);
err = abs(u_t - uexh)/abs(uexh);
assert(err < 1e-6, 'Dirichlet field error %.2e', err);

% mixed Dirichlet/Neumann, every corner mixed
nedge = numel(cgh.echnks);
isdir = mod(1:nedge, 2) == 1;

Khd = Sh + kernel.eye(0);
Khn = 0.5*kernel.eye() + kernel('helm', 'sprime', zk(1));

Kmix(nedge, nedge) = kernel();
for i = 1:nedge
    ki = Khn;
    if isdir(i), ki = Khd; end
    for j = 1:nedge
        Kmix(i,j) = ki;
    end
end

un = chnk.helm2d.kern(zk(1), srcinfo, cgh, 'sprime');
npts = arrayfun(@(c) c.npt, cgh.echnks);
dirnode = isdir(repelem(1:nedge, npts)).';
rhsmix = un;
rhsmix(dirnode) = ub(dirnode);

Amix = chunkermat(cgh, Kmix, opts);
sig = Amix \ rhsmix;
u_t = chunkerkerneval(cgh, Sh, sig, targh);
err = abs(u_t - uexh)/abs(uexh);
assert(err < 1e-6, 'mixed BC field error %.2e', err);

% 2x2 first kind system, one wavenumber per component
S2 = kernel('helm', 's', zk(2));
Zh = kernel.zeros(1,1);
Ksys(2,2) = kernel();
Ksys(1,1) = Sh; Ksys(1,2) = Zh;
Ksys(2,1) = Zh; Ksys(2,2) = S2;
Kvec = kernel(Ksys) + kernel.eye(zeros(2,2));

rhsvec = zeros(2*cgh.npt,1);
rhsvec(1:2:end) = ub;
rhsvec(2:2:end) = chnk.helm2d.kern(zk(2), srcinfo, cgh, 's');

Avec = chunkermat(cgh, Kvec, opts);
sig = Avec \ rhsvec;

Sk = {Sh, S2};
for i = 1:2
    u_t = chunkerkerneval(cgh, Sk{i}, sig(i:2:end), targh);
    uexi = chnk.helm2d.kern(zk(i), srcinfo, targinfo, 's');
    err = abs(u_t - uexi)/abs(uexi);
    assert(err < 1e-6, 'vector first kind field error %.2e (comp %d)', err, i);
end

end
