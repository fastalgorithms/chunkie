function U = wlchs_src_precomp(a,b,zsc)
%CHNK.HELSINGO.WLCHS_SRC_PRECOMP  precompute source-side data for
% wlchs_target.
%
% Inputs (panel parameterized in complex form):
%   a, b  - complex endpoints of the panel
%   zsc   - 1 x Ng complex source nodes (Gauss-Legendre on the panel)
%
% Output:
%   U     - Ng x Ng inverse Vandermonde matrix mapping function values at
%           Legendre nodes to Legendre coefficients.
%
% O(p^3) operations.  Reusable across many target points on the same panel.

cc    = (b-a)/2;
n     = numel(zsc);

% On dyadically-refined panels (RCIP) cc can shrink to ~machine eps.
% Then (zsc-mid)/cc suffers catastrophic cancellation, zsctr collapses, and
% inv(p.') is singular.  Mathematically zsctr exactly equals the real
% Gauss-Legendre nodes for a straight panel and deviates by
% O(curvature·panel_length) for a smooth curve.  When the panel is so small
% that no meaningful curvature variation is resolvable (|cc| below
% sqrt(eps)·char_scale), zsctr is unreliable; use the cached
% inverse-Vandermonde for canonical GL nodes instead.
persistent U_cache_real xs_cache n_cache
if isempty(U_cache_real) || n_cache ~= n
    [xs_gl, ~, ~, ~] = lege.exps(n);
    [p_gl, ~] = lege.pols(xs_gl(:), n-1);
    U_cache_real = inv(p_gl.');
    xs_cache = xs_gl(:);
    n_cache = n;
end

% Fall back to cache when (zsc-mid)/cc would lose precision.  Threshold:
% cc smaller than ~sqrt(eps) of the panel mid-magnitude → cancellation in
% (zsc - mid) drops too many digits.  Conservative.
ref_scale = max(1, abs((b+a)/2));
use_cache = abs(cc) < 1e-8 * ref_scale;

if use_cache
    U = U_cache_real;
else
    zsctr = (zsc(:) - (b+a)/2)/cc;
    [p,~] = lege.pols(zsctr, n-1);
    U = inv(p.');
end
end
