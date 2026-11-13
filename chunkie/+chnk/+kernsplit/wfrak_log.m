function [WfrakL, accept] = wfrak_log(ngl, trans, mscale)
%CHNK.HELSINGO.WFRAK_LOG  parameter-space log-moment correction for an
% adjacent-panel pair, used by the kernel-split close-correction
% machinery (port of testhelmos's WfrakLinit).
%
% This is the right tool when the target points lie ALONG the curve on
% an immediately-adjacent panel (i.e., on the real-axis continuation of
% the source panel in the canonical [-1,1] parameter space).  In that
% setting wlchs_target's complex-log + ifleft branch-cut machinery
% produces incorrect log moments (cuts cross through real-axis targets);
% wfrak_log uses a real-valued recurrence instead.  Use wlchs_target for
% off-curve close eval, and for the hypersingular moment HypC even for
% on-curve adjacent targets (HypC does not suffer from the same branch
% cut issue).
%
% Inputs:
%   ngl     - panel order (number of GL nodes per panel).
%   trans   - canonical-space translation: target-panel midpoint mapped
%             into source-panel parameter coordinates.  E.g. for the
%             "superdiagonal" (target panel to the LEFT of source),
%             trans = -1 - alpha; for "subdiagonal" (target to the RIGHT),
%             trans = +1 + alpha; where alpha = panel-length ratio.
%   mscale  - canonical-space scale: alpha (target-to-source size ratio).
%
% Outputs:
%   WfrakL  - na x ngl matrix.  WfrakL(i,j) is the log-moment integral
%               \int_{-1}^{1} ell_j(s) * log|s - T_i| ds
%             where ell_j is the Lagrange basis at the j-th source GL
%             node and T_i = trans + mscale*T_node_i is the i-th
%             "accepted" target's canonical position in the source frame.
%   accept  - 1 x na index vector of which target rows had |T_i|<2 (i.e.
%             which adjacent-panel target nodes are within wLCHS support).
%             Targets with |T_i|>=2 are treated as far enough that the
%             smooth GL quadrature is sufficient.

[T_nodes, ~, U_GL] = lege.exps(ngl);     % U_GL: values -> Legendre coeffs (scaled)
T_nodes = T_nodes(:);
T_targ  = trans + mscale*T_nodes;
accept  = find(abs(T_targ) < 2).';
na = numel(accept);

if na == 0
    WfrakL = zeros(0, ngl);
    return;
end

P = zeros(na, ngl+1);
Q = zeros(na, ngl);
c0 = 2*(1:ngl-1) + 1;
c1 = 1./c0;             % 1/(2k+1)
c2 = 1./(2:ngl);        % 1/(k+1)

Ta = T_targ(accept);
upp = log(abs(1 - Ta));
loo = log(abs(1 + Ta));
P(:,1) = upp - loo;
P(:,2) = 2 + Ta.*P(:,1);
for k = 1:ngl-1
    P(:,k+2) = (c0(k)*Ta.*P(:,k+1) - k*P(:,k))*c2(k);
end

Q(:,1) = 2*upp - (Ta+1).*P(:,1) - 2;
Q(:,2:ngl) = (P(:,1:ngl-1) - P(:,3:ngl+1)).*c1;

WfrakL = Q*U_GL;
end
