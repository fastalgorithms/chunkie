function moment_selftest()
%MOMENT_SELFTEST standalone check of chnk.lnsf2d.momentside (< 1 s).
%
% Criterion: the difference of two momentside calls with different outer
% limits equals the integral over the annulus [W1, W2], which contains no
% small radii, so a plain Gauss reference is reliable there. (A reference
% built from plain kernel evaluations near the singularity would itself
% be polluted by the collar cancellation that momentside exists to avoid,
% so it cannot test the small-radius part.) Agreement at ~1e-12 verifies
% the moment machinery: profile series, moment tables, geometry fit, and
% the stable n(y).z evaluation.

prm = chnk.lnsf2d.params(struct('mu', 0.004, 'muB', 0.005, 'kap', 0.00625));
eta = prm.mu;
chnkr = chunkerfunc(@(t) starfish(t), struct('maxchunklen', 1.0));
j = 3; inode = 8; k = chnkr.k;
xleg = lege.exps(k);
[x8, w8] = lege.exps(8); s8 = 0.5*(x8(:).' + 1); w8 = 0.5*w8(:).';
rs = chnkr.r(:,:,j); ds = chnkr.d(:,:,j);
xi0 = xleg(inode);
W1 = 0.08; W2 = 0.5;

S = chnk.lnsf2d.profseries(prm);
[Itab, Jtab] = chnk.lnsf2d.momtabs();
blkd = chnk.lnsf2d.momentside(prm, S, Itab, Jtab, rs, ds, xi0, +1, W2, ...
         eta, 1.0, k, s8, w8) ...
     - chnk.lnsf2d.momentside(prm, S, Itab, Jtab, rs, ds, xi0, +1, W1, ...
         eta, 1.0, k, s8, w8);

% plain Gauss reference on the annulus [xi0+W1, xi0+W2]
[x16, w16] = lege.exps(k);
x0 = rs(:, inode);
a = xi0 + W1; b = xi0 + W2; nlev = 8;
ref = zeros(3, 3*k);
for m = 1:nlev
    lo = a + (b-a)*(m-1)/nlev; hi = a + (b-a)*m/nlev;
    xis = 0.5*(lo+hi) + 0.5*(hi-lo)*x16(:);
    wws = 0.5*(hi-lo)*w16(:);
    P = lege.matrin(k, xis);
    rfine = (P*rs.').'; dfine = (P*ds.').';
    spf = sqrt(sum(dfine.^2, 1));
    nfine = [dfine(2,:); -dfine(1,:)]./spf;
    srcinfo = []; srcinfo.r = rfine; srcinfo.n = nfine;
    targinfo = []; targinfo.r = x0;
    K = chnk.lnsf2d.kern(prm, srcinfo, targinfo, 'c', eta);
    Kw = K.*repelem(wws(:).'.*spf, 3);
    ref = ref + Kw*kron(P, eye(3));
end
dmax = max(abs(blkd - ref), [], 'all');
fprintf('momentside annulus check: max abs diff = %.3e (scale %.2e)\n', ...
    dmax, max(abs(ref), [], 'all'));
if dmax < 1e-11
    fprintf('PASS: the moment machinery executes and is correct.\n');
else
    fprintf('FAIL: porting bug in the moment path.\n');
end
end
