function M = helm2d_close(type,zk,ztrg,zsrc,nz,wzp,LogC,CauC,HypC,nzt,awzp)
%CHNK.HELSINGO.HELM2D_CLOSE  close-evaluation correction matrices for
% Helmholtz layer potentials, in the kernel-split style of testhelmos.
%
% Supported types: 's' / 'single', 'd' / 'double', 'sp' / 'sprime',
% 'dp' / 'dprime' / 't'.
% These are the only types currently implemented; the LogC / CauC / HypC
% moments are kernel-independent (from wlchs_target), but the closed-form
% combination with the Helmholtz J_0 / J_1 / J_2 coefficients below is
% specific to the (i/4) H_0(kr) Green's function used by chunkie.  Adding
% another kernel requires a new case here.
%
% Output convention: M*density returns the close-evaluated layer
% potential at the targets, EXCEPT for 's' where the caller multiplies
% in awzp externally (i.e. M*(awzp.*density)).  This matches the existing
% convention used by helm2d_panel_eval.
%
% Inputs:
%   type  - 's'/'single', 'd'/'double', 'dp'/'dprime'/'t'.
%   zk    - wavenumber.
%   ztrg  - nt x 1 target points (complex form).
%   zsrc  - 1 x Ng source points on panel (complex form).
%   nz    - 1 x Ng complex unit normals at source nodes (used by 'd','dp').
%   wzp   - 1 x Ng complex contour weights z'(t)*w_GL (used by 'd').
%   LogC  - nt x Ng log-singular correction from wlchs_target (all types).
%   CauC  - nt x Ng Cauchy correction from wlchs_target ('d' only).
%   HypC  - nt x Ng hypersingular correction ('dp' only).
%   nzt   - nt x 1 complex target normals ('dp' only).
%   awzp  - 1 x Ng real arclength weights |z'(t)|*w_GL ('dp' only,
%           baked into the LogC part of M to match the 'd' convention).
%
% Output:
%   M     - nt x Ng correction matrix.

% Note: chunkie's Helmholtz S = (i/4) H_0, so the log-singular part is
% -(1/(2 pi)) log(zk r / 2) J_0(zk r) -- half of the testhelmos
% convention -- and the close-correction matrices below carry the
% corresponding factor of 1/(2 pi) instead of 1/pi.
switch lower(type)
    case {'s','single'}
        M = -besselj(0, zk*abs(ztrg - zsrc))/(2*pi) .* LogC;

    case {'d','double'}
        zdiff = ztrg - zsrc;
        tmp   = zk*abs(zdiff);
        u = tmp.*besselj(1,tmp).*imag(wzp./zdiff);
        M = -(u.*LogC + real(CauC))/(2*pi);

    case {'sp','sprime'}
        % S' = ∂_n_t G(t,s) = -(i/4) zk H_1(zk r) (n_t·(t-s))/r.
        % Helsing-Karlsson, arXiv:1711.09796, Section 4.4 (eq. KAC), with
        % chunkie's G = (i/4) H_0 introducing an extra factor 1/2 vs Helsing's
        % (i/2) H_0:
        %   G_L(r,r') = +(1/(2π)) (zk r) J_1(zk r) real(nzt/zdiff)
        %   G_C(z,τ)  = nzt · conj(nz_s) / (2π)
        %   G_H       = 0
        % The Cauchy moment G_C is absorbed into G_0 on smooth Γ (Helsing
        % demo12.m drops it), but is required near corners and branch points.
        % The wcmpC delta from wlchs_target gives the per-pair (analytical -
        % smooth GL) Cauchy correction.
        % besselj(1,·) form keeps the formula valid for real OR complex zk
        % (Helsing's imag(K) trick assumes real zk).
        if nargin < 8 || isempty(CauC)
            error('CHNK.HELSINGO.HELM2D_CLOSE: ''sp'' requires CauC');
        end
        if nargin < 10 || isempty(nzt)
            error('CHNK.HELSINGO.HELM2D_CLOSE: ''sp'' requires target normals nzt');
        end
        if nargin < 11 || isempty(awzp)
            error('CHNK.HELSINGO.HELM2D_CLOSE: ''sp'' requires source weights awzp');
        end
        nzt  = nzt(:);             % nt x 1
        nz   = nz(:).';            % 1 x Ng
        awzp = awzp(:).';          % 1 x Ng
        zdiff = ztrg - zsrc;        % nt x Ng
        tmp = zk*abs(zdiff);
        cosNT = real(nzt ./ zdiff);
        Mlog = (tmp .* besselj(1,tmp) .* cosNT) .* LogC .* awzp / (2*pi);
        Mcau = real(nzt .* conj(nz) .* CauC) / (2*pi);
        M = Mlog + Mcau;

    case {'dp','dprime','t'}
        % T = D' (hypersingular, normal-derivative of double-layer).
        % Translated from testhelmos ToperC_near, divided by 2 (chunkie
        % convention). awzp is baked into the LogC piece so that
        % M*density is the close correction.
        if nargin < 9 || isempty(HypC)
            error('CHNK.HELSINGO.HELM2D_CLOSE: ''dp'' requires HypC');
        end
        if nargin < 10 || isempty(nzt)
            error('CHNK.HELSINGO.HELM2D_CLOSE: ''dp'' requires target normals nzt');
        end
        if nargin < 11 || isempty(awzp)
            error('CHNK.HELSINGO.HELM2D_CLOSE: ''dp'' requires source weights awzp');
        end
        nzt  = nzt(:);            % nt x 1
        nz   = nz(:).';           % 1 x Ng
        awzp = awzp(:).';         % 1 x Ng
        zdiff = ztrg - zsrc;      % nt x Ng (broadcast)
        tmp = zk*abs(zdiff);
        % Coefficient of LogC: smooth function of (target,source).
        Ta1 = (zk^2) * besselj(1,tmp)./tmp .* real(nzt .* conj(nz));
        D1  = -real(nz   ./ zdiff);
        D2  =  real(nzt  ./ zdiff);
        Tb1 = tmp.^2 .* besselj(2,tmp) .* D1 .* D2;
        Mlog = -(Ta1 + Tb1) .* LogC / (2*pi) .* awzp;
        Mhyp = -real(nzt .* HypC) / (2*pi);
        M = Mlog + Mhyp;

    otherwise
        error('CHNK.HELSINGO.HELM2D_CLOSE: unsupported type %s', type);
end
end
