chunkermat_lnsf2dNearquadTest0();

function chunkermat_lnsf2dNearquadTest0()
%CHUNKERMAT_LNSF2DNEARQUADTEST
%
% Tests for the window-quadrature correction (chnk.lnsf2d.nearcorrect) of
% self/near blocks of the thermoviscous combined-layer operator on chunks
% that are long compared to the boundary-layer decay length:
%
%   (1) block-level unit test: corrected self blocks agree with a
%       brute-force dyadically-refined reference (with stable difference-
%       form displacements) to ~1e-10 on a coarse discretization
%   (2) MMS comparison, plain chunkermat vs corrected, as a function of
%       the chunk length in units of delta = 1/Im(kh); also checks that
%       the correction is a no-op when no chunk is flagged
%
% see also CHNK.LNSF2D.NEARCORRECT, CHNK.LNSF2D.FLAGCHUNKS

iseed = 8675309;
rng(iseed);

prm = chnk.lnsf2d.params();
eta = 1.0;
delta = 1/imag(prm.kh);
kabs = max(abs(prm.kh), abs(prm.kv));
fprintf('delta = %.3f, |kmax| = %.3f\n', delta, kabs);

fkern = kernel(@(s,t) chnk.lnsf2d.kern(prm, s, t, 'c', eta));
fkern.opdims = [3 3];
fkern.sing = 'log';

narms = 3; amp = 0.25;

% ---------------------------------------------------------------
% (1) block-level unit test on a coarse discretization
% ---------------------------------------------------------------
nch = 6;
chnkr = chunkerfuncuni(@(t) starfish(t, narms, amp), nch);
flags = chnk.lnsf2d.flagchunks(chnkr, prm);
fprintf('nch = %d: flagged chunks: %s\n', nch, mat2str(flags));
assert(~isempty(flags));

sysmat = chunkermat(chnkr, fkern);
syscorr = chnk.lnsf2d.nearcorrect(sysmat, chnkr, prm, eta);

k = chnkr.k;
kfun = @(s, t) chnk.lnsf2d.kern(prm, s, t, 'c', eta);
maxdiff_ref = 0; maxdiff_ggq = 0;
for j = flags(1:min(2, numel(flags)))
    for inode = [1, 6, 16]           % includes both edge nodes
        i = (j-1)*k + inode;
        irow = 3*(i-1) + (1:3);
        icol = 3*k*(j-1) + (1:3*k);
        blkref = brute_self_block(chnkr, prm, kfun, j, inode);
        e1 = max(abs(syscorr(irow, icol) - blkref), [], 'all');
        e2 = max(abs(sysmat(irow, icol) - blkref), [], 'all');
        maxdiff_ref = max(maxdiff_ref, e1);
        maxdiff_ggq = max(maxdiff_ggq, e2);
    end
end
fprintf('corrected vs brute reference: %.3e   (plain ggq vs ref: %.3e)\n', ...
    maxdiff_ref, maxdiff_ggq);
assert(maxdiff_ref < 1e-9);

% ---------------------------------------------------------------
% (2) MMS comparison, plain vs corrected
% ---------------------------------------------------------------
srcpts = []; srcpts.r = [0.10, -0.20; -0.35, 0.15];
strengths = randn(6,1) + 1i*randn(6,1);
tt = 2*pi*(0:7)/8;
targs = [2.2*cos(tt), 1.62*cos(tt(1:2:end)+0.3); ...
         2.2*sin(tt), 1.62*sin(tt(1:2:end)+0.3)];
uex = @(pts) chnk.lnsf2d.kern(prm, srcpts, struct('r', pts), 's')*strengths;
utrue = uex(targs);

% correction options: flag more aggressively and always recompute the
% full adjacent (near-rule) blocks of flagged chunks
copts = []; copts.tau = 2.0;

nchs = [24, 12, 8, 6, 4];
res = zeros(2, length(nchs));
ifar = 1:8; inear8 = 9:12;    % target ring indices (8 far + 4 near)
for imc = 1:length(nchs)
    chnkr = chunkerfuncuni(@(t) starfish(t, narms, amp), nchs(imc));
    maxlen = max(chunklen(chnkr));
    flags = chnk.lnsf2d.flagchunks(chnkr, prm, copts.tau);
    sysmat = chunkermat(chnkr, fkern);
    syscorr = chnk.lnsf2d.nearcorrect(sysmat, chnkr, prm, eta, copts);
    if isempty(flags)
        assert(norm(syscorr - sysmat, 'fro') == 0);   % no-op when unflagged
    end
    bdrypts = reshape(chnkr.r, 2, chnkr.k*chnkr.nch);
    rhs = uex(bdrypts);
    opts = []; opts.usesmooth = false;
    psis = cell(2,1);
    for iv = 1:2
        if iv == 1, sys = 0.5*eye(size(sysmat,1)) + sysmat;
        else, sys = 0.5*eye(size(sysmat,1)) + syscorr; end
        psis{iv} = sys \ rhs;
        usol = chunkerkerneval(chnkr, fkern, psis{iv}, targs, opts);
        res(iv, imc) = max(abs(usol - utrue))/max(abs(utrue));
        if iv == 2
            % diagnostics: error split by target ring, plus an independent
            % evaluator (graded-piece quadrature, kernel-agnostic)
            errt = abs(usol - utrue)/max(abs(utrue));
            errt = reshape(errt, 3, []);
            efar = max(errt(:, ifar), [], 'all');
            enear = max(errt(:, inear8), [], 'all');
            usol2 = evalref(chnkr, prm, eta, psis{iv}, targs);
            e2 = max(abs(usol2 - utrue))/max(abs(utrue));
            fprintf(['      corrected: far-ring %.3e  near-ring %.3e  ' ...
                'ref-eval %.3e  |psi_c - psi_p| = %.3e\n'], ...
                efar, enear, e2, max(abs(psis{2} - psis{1})));
        end
    end
    fprintf('nch = %3d  maxlen/delta = %5.2f  nflag = %2d  plain = %.3e  corrected = %.3e\n', ...
        nchs(imc), maxlen/delta, numel(flags), res(1, imc), res(2, imc));
end
% corrected matrix should not be worse than plain
assert(all(res(2,:) < 2*res(1,:) + 1e-10));

% ---------------------------------------------------------------
% (3) decouple boundary-layer scale from geometry resolution:
% nu, kappa decreased 10x (delta -> delta/sqrt(10), |kv|/|ka| ~ 16).
% With the corrected quadrature the error should be set by the number
% of chunks needed for the geometry/acoustic scales, NOT by
% maxchunklen/delta: chunks 8-10 delta long should stay at ~1e-10.
% ---------------------------------------------------------------
popts = []; popts.mu = 0.004; popts.muB = 0.005; popts.kap = 0.00625;
prm10 = chnk.lnsf2d.params(popts);
delta10 = 1/imag(prm10.kh);
fprintf('\nnu/10 sweep: delta = %.4f, |kv|/|ka| = %.1f\n', ...
    delta10, abs(prm10.kv)/abs(prm10.ka));
fkern10 = kernel(@(s,t) chnk.lnsf2d.kern(prm10, s, t, 'c', eta));
fkern10.opdims = [3 3];
fkern10.sing = 'log';
uex10 = @(pts) chnk.lnsf2d.kern(prm10, srcpts, struct('r', pts), 's')*strengths;
utrue10 = uex10(targs);
nchs10 = [32, 24, 16, 12, 10];
res10 = zeros(2, length(nchs10));
for imc = 1:length(nchs10)
    chnkr = chunkerfuncuni(@(t) starfish(t, narms, amp), nchs10(imc));
    maxlen = max(chunklen(chnkr));
    sysmat = chunkermat(chnkr, fkern10);
    syscorr = chnk.lnsf2d.nearcorrect(sysmat, chnkr, prm10, eta, copts);
    bdrypts = reshape(chnkr.r, 2, chnkr.k*chnkr.nch);
    rhs = uex10(bdrypts);
    opts = []; opts.usesmooth = false;
    for iv = 1:2
        if iv == 1, sys = 0.5*eye(size(sysmat,1)) + sysmat;
        else, sys = 0.5*eye(size(sysmat,1)) + syscorr; end
        psi = sys \ rhs;
        usol = chunkerkerneval(chnkr, fkern10, psi, targs, opts);
        res10(iv, imc) = max(abs(usol - utrue10))/max(abs(utrue10));
    end
    fprintf('nch = %3d  maxlen/delta = %5.2f  plain = %.3e  corrected = %.3e\n', ...
        nchs10(imc), maxlen/delta10, res10(1, imc), res10(2, imc));
end
assert(res10(2, nchs10 == 12) < 1e-8);

% ---------------------------------------------------------------
% (4) realistic viscosity: mu = 4e-6, kappa = 6.25e-6, so that
% |kv|/|ka| = 500 and delta ~ 3e-3 (chunks are 170-220 delta long).
% The boundary layer is far below the grading scale of the ggq aux
% nodes, so the plain rules fail at ~1e-4 while the corrected blocks
% stay at the geometry/data-limited floor (~1e-10). The CFIE coupling
% must be rescaled, eta = mu, since the single-layer blocks scale like
% 1/mu (mobility) and 1/kappa.
% ---------------------------------------------------------------
popts = []; popts.mu = 4e-6; popts.muB = 5e-6; popts.kap = 6.25e-6;
prmr = chnk.lnsf2d.params(popts);
deltar = 1/imag(prmr.kh);
etar = prmr.mu;
fprintf('\nrealistic-nu sweep: delta = %.2e, |kv|/|ka| = %.0f, eta = mu\n', ...
    deltar, abs(prmr.kv)/abs(prmr.ka));
fkernr = kernel(@(s,t) chnk.lnsf2d.kern(prmr, s, t, 'c', etar));
fkernr.opdims = [3 3];
fkernr.sing = 'log';
uexr = @(pts) chnk.lnsf2d.kern(prmr, srcpts, struct('r', pts), 's')*strengths;
utruer = uexr(targs);
nchsr = [16, 12];
resr = zeros(2, length(nchsr));
for imc = 1:length(nchsr)
    chnkr = chunkerfuncuni(@(t) starfish(t, narms, amp), nchsr(imc));
    maxlen = max(chunklen(chnkr));
    sysmat = chunkermat(chnkr, fkernr);
    syscorr = chnk.lnsf2d.nearcorrect(sysmat, chnkr, prmr, etar, copts);
    bdrypts = reshape(chnkr.r, 2, chnkr.k*chnkr.nch);
    rhs = uexr(bdrypts);
    opts = []; opts.usesmooth = false;
    for iv = 1:2
        if iv == 1, sys = 0.5*eye(size(sysmat,1)) + sysmat;
        else, sys = 0.5*eye(size(sysmat,1)) + syscorr; end
        psi = sys \ rhs;
        usol = chunkerkerneval(chnkr, fkernr, psi, targs, opts);
        resr(iv, imc) = max(abs(usol - utruer))/max(abs(utruer));
    end
    fprintf('nch = %3d  maxlen/delta = %6.1f  plain = %.3e  corrected = %.3e\n', ...
        nchsr(imc), maxlen/deltar, resr(1, imc), resr(2, imc));
end
assert(all(resr(2,:) < 1e-8));
assert(resr(1,1)/resr(2,1) > 1e3);   % correction should win by >3 digits

end

% ------------------------------------------------------------------
function u = evalref(chnkr, prm, eta, psi, targs)
%EVALREF independent representation evaluator: per (target, chunk),
% graded pieces with |kmax| x arclength <= 3 in the decay window and 2:1
% grading toward the closest point, plain 16-pt Gauss per piece.
k = chnkr.k;
nch = chnkr.nch;
xleg = lege.exps(k);
[x16, w16] = lege.exps(k);
kfun = @(s, t) chnk.lnsf2d.kern(prm, s, t, 'c', eta);
kabs = max(abs(prm.kh), abs(prm.kv));
kim = min(imag(prm.kh), imag(prm.kv));
nt = size(targs, 2);
u = zeros(3*nt, 1);
for it = 1:nt
    x0 = targs(:, it);
    val = zeros(3, 1);
    for j = 1:nch
        rs = chnkr.r(:,:,j); ds = chnkr.d(:,:,j);
        sp = sqrt(sum(ds.^2, 1));
        spmax = max(sp); spmin = min(sp);
        hbl = 3.0/(kabs*spmax);
        W = 34.0/(kim*spmin);
        [dmin, imin] = min((rs(1,:)-x0(1)).^2 + (rs(2,:)-x0(2)).^2);
        dmin = sqrt(dmin);
        dfloor = max(dmin/(2*spmax), 1e-12);
        % graded pieces toward the closest node parameter
        xic = xleg(imin);
        segs = zeros(0, 2); stack = [-1.0, 1.0];
        while ~isempty(stack)
            lo = stack(end,1); hi = stack(end,2); stack(end,:) = [];
            ln = hi - lo;
            if xic < lo, dist = lo - xic;
            elseif xic > hi, dist = xic - hi;
            else, dist = 0.0; end
            graded = ln <= 2.0*dist + 1e-15;
            ok = graded && ((ln <= hbl) || (dist >= W && ln <= 0.5*dist + 1e-15));
            if ok || ln <= dfloor
                segs(end+1, :) = [lo, hi]; %#ok<AGROW>
            else
                mid = 0.5*(lo + hi);
                stack(end+1, :) = [lo, mid]; %#ok<AGROW>
                stack(end+1, :) = [mid, hi]; %#ok<AGROW>
            end
        end
        xis = []; wws = [];
        for m = 1:size(segs, 1)
            lo = segs(m,1); hi = segs(m,2);
            xis = [xis; 0.5*(lo+hi) + 0.5*(hi-lo)*x16(:)]; %#ok<AGROW>
            wws = [wws; 0.5*(hi-lo)*w16(:)]; %#ok<AGROW>
        end
        P = lege.matrin(k, xis);
        rfine = (P*rs.').';
        dfine = (P*ds.').';
        spf = sqrt(sum(dfine.^2, 1));
        nfine = [dfine(2,:); -dfine(1,:)]./spf;
        srcinfo = []; srcinfo.r = rfine; srcinfo.n = nfine;
        targinfo = []; targinfo.r = x0;
        K = kfun(srcinfo, targinfo);
        Kw = K.*repelem(wws(:).'.*spf, 3);
        psij = psi(3*k*(j-1) + (1:3*k));
        val = val + Kw*(kron(P, eye(3))*psij);
    end
    u(3*(it-1)+(1:3)) = val;
end
end

function blk = brute_self_block(chnkr, prm, kfun, j, inode)
%BRUTE_SELF_BLOCK reference (3 x 3k) self block for target node inode on
% chunk j: dyadic refinement toward the target in the chunk parameter,
% with displacements computed in stable difference form through the
% chunk's derivative polynomial.
k = chnkr.k;
rs = chnkr.r(:,:,j); ds = chnkr.d(:,:,j);
[xleg, wleg] = lege.exps(k);
[x8, w8] = lege.exps(8);
s8 = 0.5*(x8(:).' + 1); w8 = 0.5*w8(:).';
xi0 = xleg(inode);
% dyadic partition of [-1,1] toward xi0
segs = zeros(0,2);
for side = 1:2
    if side == 1, lo = -1.0; hi = xi0; tolo = false;
    else, lo = xi0; hi = 1.0; tolo = true; end
    len = hi - lo; d = 0;
    while d < 46 && len > 2e-14
        len = len/2;
        if tolo, segs(end+1,:) = [lo+len, lo+2*len]; %#ok<AGROW>
        else, segs(end+1,:) = [hi-2*len, hi-len]; end %#ok<AGROW>
        d = d+1;
    end
    if tolo, segs(end+1,:) = [lo, lo+len]; %#ok<AGROW>
    else, segs(end+1,:) = [hi-len, hi]; end %#ok<AGROW>
end
xis = []; wws = [];
for m = 1:size(segs,1)
    lo = segs(m,1); hi = segs(m,2);
    xis = [xis; 0.5*(lo+hi) + 0.5*(hi-lo)*xleg(:)]; %#ok<AGROW>
    wws = [wws; 0.5*(hi-lo)*wleg(:)]; %#ok<AGROW>
end
keep = abs(xis - xi0) > 1e-14;
xis = xis(keep); wws = wws(keep);
nq = numel(xis);
P = lege.matrin(k, xis);
dfine = (P*ds.').';
spf = sqrt(sum(dfine.^2, 1));
nfine = [dfine(2,:); -dfine(1,:)]./spf;
dxi = xis(:).' - xi0;
taus = xi0 + dxi(:)*s8;
P8 = lege.matrin(k, taus(:));
dtau = reshape((P8*ds.').', 2, nq, 8);
zint = zeros(2, nq);
for m = 1:numel(w8)
    zint = zint + w8(m)*dtau(:,:,m);
end
z = -dxi.*zint;
srcinfo = []; srcinfo.r = -z; srcinfo.n = nfine;
targinfo = []; targinfo.r = [0; 0];
K = kfun(srcinfo, targinfo);
Kw = K.*repelem(wws(:).'.*spf, 3);
blk = Kw*kron(P, eye(3));
end
