function sysmat = nearcorrect(sysmat, chnkr, prm, eta, opts)
%CHNK.LNSF2D.NEARCORRECT recompute self/near matrix blocks of the
% thermoviscous combined-layer operator for chunks that are long compared
% to the boundary-layer decay length, where the standard GGQ rules lose
% accuracy.
%
% For every flagged source chunk (see CHNK.LNSF2D.FLAGCHUNKS) and every
% target node on or near it, the corresponding (3 x 3k) block is
% recomputed with a window quadrature whose node count is independent of
% kh, kv:
%
%   * a small "singular piece" around the target, sized so that
%     |k_max| h0 ~ H0FAC, integrated by moment integration
%     (CHNK.LNSF2D.MOMENTSIDE): the kernel's radial log-series is
%     contracted against exact moments of a polynomial fit of the
%     geometry, so the kernel's large collar values (up to 1/delta^2 in
%     the temperature rows) are never sampled;
%   * plain Gauss pieces with |k_max| x arclength <= HFAC covering the
%     boundary-layer decay window (Im k x dist <= WFAC), with 2:1 grading
%     toward the target for the log tail;
%   * 2:1 graded pieces outside the window.
%
% Syntax:
%   sysmat = chnk.lnsf2d.nearcorrect(sysmat, chnkr, prm)
%   sysmat = chnk.lnsf2d.nearcorrect(sysmat, chnkr, prm, eta, opts)
%
% Input:
%   sysmat - system matrix from chunkermat(chnkr, kernel('c'-type lnsf2d))
%   chnkr - chunker object
%   prm - parameter struct from chnk.lnsf2d.params
%   eta - combined-field coupling (default 1.0); must match the kernel
%         used to build sysmat
%   opts - options struct:
%       opts.tau (3.0) flagging threshold, see flagchunks
%       opts.h0fac (4.0), opts.hfac (3.0), opts.wfac (34.0) rule shape
%       opts.reach (1.0) targets within reach*chunklen are recomputed
%       opts.scale (1.0) overall kernel scaling (2 for the rcip
%           convention I + 2(D + eta S))
%       opts.tscale (1.0) thermal-density conjugation; must match the
%           kernel used to build sysmat
%       opts.skipsets ({}) cell array of node-index sets (one per rcip
%           vertex, e.g. unique(ceil(rcipsav{v}.starind/3))): a pair is
%           skipped iff the target node and all nodes of the source chunk
%           lie in the same set (those blocks are owned by inv(R) - I)
%
% Output:
%   sysmat - corrected system matrix
%
% see also CHNK.LNSF2D.KERN, CHNK.LNSF2D.FLAGCHUNKS, CHNK.LNSF2D.MOMENTSIDE

if nargin < 4 || isempty(eta), eta = 1.0; end
if nargin < 5, opts = []; end
tau = getf(opts, 'tau', 3.0);
h0fac = getf(opts, 'h0fac', 4.0);
hfac = getf(opts, 'hfac', 3.0);
wfac = getf(opts, 'wfac', 34.0);
reach = getf(opts, 'reach', 1.0);
scale = getf(opts, 'scale', 1.0);
skipsets = getf(opts, 'skipsets', {});
frmin = 0.12;

k = chnkr.k;
nch = chnkr.nch;
xleg = lege.exps(k);
[x8, w8] = lege.exps(8);
s8 = 0.5*(x8(:).' + 1); w8 = 0.5*w8(:).';

kabs = max(abs(prm.kh), abs(prm.kv));
kim = min(imag(prm.kh), imag(prm.kv));

flags = chnk.lnsf2d.flagchunks(chnkr, prm, tau);
if isempty(flags), return; end

lens = chunklen(chnkr);
allr = reshape(chnkr.r, 2, k*nch);
tscale = getf(opts, 'tscale', 1.0);   % must match the kernel used to
                                      % build the system matrix
kfun = @(s, t) chnk.lnsf2d.kern(prm, s, t, 'c', eta, tscale);
Sprof = chnk.lnsf2d.profseries(prm);   % radial profile series
[Itab, Jtab] = chnk.lnsf2d.momtabs(); % exact moment tables

for j = flags
    rs = chnkr.r(:,:,j); ds = chnkr.d(:,:,j);
    sp = sqrt(sum(ds.^2, 1));
    spmax = max(sp); spmin = min(sp);
    hbl = hfac/(kabs*spmax);
    W = wfac/(kim*spmin);
    % targets: nodes of chunk j, ALL nodes of adjacent chunks (these are
    % the pairs chunkermat treated with the ggq near rule, whose accuracy
    % also degrades on long chunks), plus nodes of other chunks within
    % reach of chunk j
    dd1 = allr(1,:).' - rs(1,:); dd2 = allr(2,:).' - rs(2,:);
    mindist = sqrt(min(dd1.^2 + dd2.^2, [], 2)).';
    inear = find(mindist < reach*lens(j));
    adjch = chnkr.adj(:, j); adjch = adjch(adjch > 0);
    for a = adjch(:).'
        inear = union(inear, (a-1)*k + (1:k));
    end
    inear = inear(:).';
    % exact rcip-star skip: is the source chunk contained in a star set?
    jstars = [];
    if ~isempty(skipsets)
        jnodes = (j-1)*k + (1:k);
        for v = 1:numel(skipsets)
            if ~isempty(skipsets{v}) && all(ismember(jnodes, skipsets{v}))
                jstars(end+1) = v; %#ok<AGROW>
            end
        end
    end
    for i = inear
        ich = floor((i-1)/k) + 1;
        inode = i - (ich-1)*k;
        x0 = allr(:, i);
        if ~isempty(jstars)
            isskip = false;
            for v = jstars
                if ismember(i, skipsets{v}), isskip = true; break; end
            end
            if isskip, continue; end
        end
        if ich == j
            % ----- self target: moment-integrated singular piece + bulk -----
            xi0 = xleg(inode);
            sp0 = sp(inode);
            h0 = h0fac/(kabs*sp0);
            [alpha, h0] = singpiece(xi0, -1.0, 1.0, h0, xleg, frmin);
            beta = alpha + h0;
            % the kernel's radial log-series is contracted against exact
            % moments of a polynomial fit of the geometry, so the kernel's
            % large collar values (up to 1/delta^2 in the T rows) are
            % never sampled; sampling them costs eps times that amplitude
            % no matter what rule is used
            blk = chnk.lnsf2d.momentside(prm, Sprof, Itab, Jtab, ...
                    rs, ds, xi0, -1.0, xi0 - alpha, eta, tscale, k, s8, w8) ...
                + chnk.lnsf2d.momentside(prm, Sprof, Itab, Jtab, ...
                    rs, ds, xi0, +1.0, beta - xi0, eta, tscale, k, s8, w8);
            % bulk + tail pieces
            segs = [bulkpieces(-1.0, alpha, xi0, hbl, W, 0.02*h0); ...
                    bulkpieces(beta, 1.0, xi0, hbl, W, 0.02*h0)];
            blk = blk + segblock(kfun, rs, ds, xleg, segs, x0, k);
        else
            % ----- near target on another chunk: graded bulk only -----
            [dmin, imin] = min((rs(1,:)-x0(1)).^2 + (rs(2,:)-x0(2)).^2);
            dmin = sqrt(dmin);
            xic = xleg(imin);
            dfloor = max(dmin/(2*spmax), 1e-12);
            segs = bulkpieces(-1.0, 1.0, xic, hbl, W, dfloor);
            blk = segblock(kfun, rs, ds, xleg, segs, x0, k);
        end
        irow = 3*(i-1) + (1:3);
        icol = 3*k*(j-1) + (1:3*k);
        sysmat(irow, icol) = scale*blk;
    end
end

end

% ------------------------------------------------------------------
function blk = segblock(kfun, rs, ds, xleg, segs, x0, k)
% plain Gauss on segments segs (nseg x 2) in xi, target at x0
blk = zeros(3, 3*k);
if isempty(segs), return; end
[x16, w16] = lege.exps(k);
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
blk = Kw*kron(P, eye(3));
end

function [alpha, h0] = singpiece(xi0, a, b, h0, xleg, frmin)
% choose [alpha, alpha+h0] in [a,b] with xi0 at Legendre node j0, the
% piece extending at least frmin*h0 on each side of xi0
h0 = min(h0, b - a);
fr = (1.0 + xleg(:).')/2;
interior = (fr >= frmin) & (fr <= 1 - frmin);
lo_need = xi0 - a; hi_need = b - xi0;
feas = interior & (fr*h0 <= lo_need) & ((1-fr)*h0 <= hi_need);
if any(feas)
    cand = find(feas);
    [~, m] = min(abs(fr(cand) - 0.5));
    j0 = cand(m);
else
    cand = find(interior);
    h0s = min(lo_need./fr(cand), hi_need./(1 - fr(cand)));
    [h0m, m] = max(h0s);
    j0 = cand(m);
    h0 = 0.95*min(h0m, h0);
end
alpha = xi0 - fr(j0)*h0;
end

function segs = bulkpieces(a, b, xic, hbl, W, dfloor)
% cover [a,b] with pieces: len <= 2*dist(xic) always (log-tail grading),
% len <= hbl inside the decay window, len <= max(dist/2, hbl) outside;
% never split below dfloor. Iterative 2:1 bisection.
segs = zeros(0, 2);
if b <= a + 1e-15, return; end
stack = [a, b];
while ~isempty(stack)
    lo = stack(end, 1); hi = stack(end, 2);
    stack(end, :) = [];
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
end

function v = getf(o, f, dflt)
if isfield(o, f), v = o.(f); else, v = dflt; end
end
