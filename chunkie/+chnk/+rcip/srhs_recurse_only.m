function r_f_star = srhs_recurse_only(rcipsav_corner, srhs_eval)
%CHNK.RCIP.SRHS_RECURSE_ONLY  Re-run the singular-RHS forward recursion at a
% single corner using cached level data, WITHOUT rebuilding any system matrix.
%
% Use this when the BIE operator is fixed (cached in rcipsav by a prior
% chunkermat call with opts.cache_mat_for_srhs=true) but the right-hand
% side data_in changes — e.g., when computing an I2I matrix by solving the
% same BIE for many input modes. Speeds up by skipping the per-iter
% rebuild of bulk + corner-level matrices.
%
% Inputs:
%   rcipsav_corner - per-corner rcipsav cell entry produced by Rcompchunk
%                    with cache_mat_for_srhs=true. Required fields:
%                      .nsub, .nedge, .ndim, .Pbc, .PWbc,
%                      .starL, .circL, .starS, .circS, .use_myinv,
%                      .R{1..nsub+1}, .MAT_full{1..nsub},
%                      .chnkrlocals_array{1..nsub}, .ctr,
%                      .vert0, .iedgechunks
%   srhs_eval      - function handle (same signature as in Rcompchunk):
%                      b_ib = srhs_eval(chnkrlocal_global, vert0,
%                                       edge_indices, level, ndim)
%
% Output:
%   r_f_star       - 2*k*ndim*nedge vector at the corner.
%
% See chnk.rcip.Rfstep, chnk.rcip.Rcompchunk.

nsub      = rcipsav_corner.nsub;
nedge     = rcipsav_corner.nedge;
ndim      = rcipsav_corner.ndim;
Pbc       = rcipsav_corner.Pbc;
PWbc      = rcipsav_corner.PWbc;
starL     = rcipsav_corner.starL;
circL     = rcipsav_corner.circL;
starS     = rcipsav_corner.starS;
circS     = rcipsav_corner.circS;
use_myinv = rcipsav_corner.use_myinv;
ctr       = rcipsav_corner.ctr;
vert0     = rcipsav_corner.vert0;
iedgechunks = rcipsav_corner.iedgechunks;

r_f_star = [];
for level = 1:nsub
    R_lev   = rcipsav_corner.R{level};
    MAT_lev = rcipsav_corner.MAT_full{level};
    chnkrlocal_lev = rcipsav_corner.chnkrlocals_array{level};

    % Translate local mesh to global coords for srhs_eval.
    chnkrlocal_global = chnkrlocal_lev;
    for ie = 1:nedge
        chnkrlocal_global(ie) = chnkrlocal_lev(ie) + ctr(:,ie);
    end
    b_ib = srhs_eval(chnkrlocal_global, vert0, iedgechunks(1,:), level, ndim);
    b_ib = b_ib(:);

    if level == 1
        r_f_star = R_lev * b_ib(starL);
    end
    r_f_star = chnk.rcip.Rfstep(MAT_lev, R_lev, r_f_star, b_ib(circL), ...
                                Pbc, PWbc, starL, circL, starS, circS, use_myinv);
end
end
