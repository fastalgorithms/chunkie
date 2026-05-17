function rhs_mod = srhs_modify_rhs(rhs_bare, rcipsav)
%CHNK.RCIP.SRHS_MODIFY_RHS apply singular-RHS correction at corner-star.
%
% In the singular-RHS RCIP framework (Helsing & Karlsson 2022), expressed
% in chunkie's hat_rho convention where sol = R*v_tilde + r_f_star, the
% only modification to the standard amat\rhs pipeline is at corner-star
% indices: the bare RHS values are replaced by inv(R) * r_f_star_corner.
%
% Mathematically (paper eq. for disc2, rewritten for hat_rho):
%   amat * hat_rho = rhs
%     amat(cs, cs)    = inv(R)
%     amat(cs, ~cs)   = K(cs, ~cs)
%     amat(~cs, ~cs)  = I + K(~cs, ~cs)
%   rhs(~cs) = b_coarse(~cs)               (unchanged)
%   rhs( cs) = inv(R) * r_f_star_corner    (replaced)
%
% The far-from-corner entries of b_coarse are correct as-is from a bare
% chunkermat(rhs_kerns) * data_in evaluation; only the corner-star entries
% need the R_f-recursion result.
%
% Inputs:
%   rhs_bare - bare RHS vector (e.g. chunkermat(rhs_kerns,opts_no_rcip) * data_in).
%   rcipsav  - per-corner cell array returned as the third output of
%              chunkermat when opts.srhs_eval is set.
%
% Output:
%   rhs_mod - RHS vector with corner-star entries replaced.
%
% See also: chnk.rcip.Rcompchunk, chnk.rcip.Rfstep, chunkermat.

rhs_mod = rhs_bare;
for ivert = 1:numel(rcipsav)
    rcs = rcipsav{ivert};
    if isempty(rcs); continue; end
    if ~isfield(rcs, 'r_f_star') || isempty(rcs.r_f_star); continue; end
    if ~isfield(rcs, 'starind') || isempty(rcs.starind); continue; end

    nsub = rcs.nsub;
    R_final = rcs.R{nsub+1};
    rhs_mod(rcs.starind) = R_final \ rcs.r_f_star;
end
end
