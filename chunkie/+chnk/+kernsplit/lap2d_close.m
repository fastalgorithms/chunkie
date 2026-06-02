function M = lap2d_close(type, ztrg, zsrc, nz, wzp, LogC, CauC, HypC, nzt, awzp)
%CHNK.HELSINGO.LAP2D_CLOSE  close-evaluation correction matrices for 2D
% Laplace layer potentials, in the kernel-split style of Helsing & Ojala.
%
% Conventions (chunkie's): G = -log|r|/(2*pi); n_y is outward-pointing
% complex normal at source; |dy| = arclength element.
%
% Layer types implemented:
%   's'/'single'   - SLP   S(tau)(x) = -1/(2pi) integral log|x-y| tau ds
%   'd'/'double'   - DLP   D(tau)(x) =  1/(2pi) integral (n_y . (x-y))/r^2 tau ds
%   'sp'/'sprime'  - normal derivative of SLP at target
%   'dp'/'dprime'  - normal derivative of DLP at target (hypersingular)
%
% Output convention follows helm2d_close: M*density gives the close-
% evaluated layer potential at targets, EXCEPT for 's' where the caller
% multiplies in awzp externally (i.e. M*(awzp.*density)).
%
% Inputs:
%   type  - 's'/'single', 'd'/'double', 'sp'/'sprime', 'dp'/'dprime'.
%   ztrg  - nt x 1 complex target points.
%   zsrc  - 1 x Ng complex source nodes.
%   nz    - 1 x Ng complex unit normals at sources (for 'd', 'dp').
%   wzp   - 1 x Ng complex contour weights z'(t)*w_GL (for sp/dp).
%   LogC  - nt x Ng log-singular correction from wlchs_target.
%   CauC  - nt x Ng Cauchy correction from wlchs_target ('d', 'sp').
%   HypC  - nt x Ng hypersingular correction ('dp').
%   nzt   - nt x 1 complex target normals ('sp', 'dp').
%   awzp  - 1 x Ng real arclength weights ('dp', baked into LogC piece).
%
% Output:
%   M     - nt x Ng correction matrix.

switch lower(type)
    case {'s','single'}
        % S = -log|r|/(2pi).  Per derivation: M = -LogC/(2pi) and caller
        % multiplies M*(awzp.*density) externally.
        M = -LogC/(2*pi);

    case {'d','double'}
        % D = (1/2pi) Re(n_y/(z-w)).  Cauchy moment 1i*CauC corresponds to
        % I_C correction; its real part divided by 2pi (with sign) is the
        % matrix correction.  Verified by zk -> 0 limit of helm2d_close.
        if nargin < 7 || isempty(CauC)
            error('CHNK.HELSINGO.LAP2D_CLOSE: ''d'' requires CauC');
        end
        M = -real(CauC)/(2*pi);

    case {'sp','sprime'}
        % S' = n_x . grad_x S = (1/2pi) Re(n_x * I_C[tau/(in_y)]).  The
        % wcmpC delta moment for source g = tau/(in_y) is wcmpC[i,j]/nz_j,
        % so the correction is:
        %   M_correction = (1/2pi) Re( n_x_i * wcmpC[i,j] / nz_j )
        if nargin < 7 || isempty(CauC)
            error('CHNK.HELSINGO.LAP2D_CLOSE: ''sp'' requires CauC');
        end
        if nargin < 9 || isempty(nzt)
            error('CHNK.HELSINGO.LAP2D_CLOSE: ''sp'' requires target normals nzt');
        end
        nzt = nzt(:);
        nz  = nz(:).';
        M = real(nzt .* (CauC ./ nz)) / (2*pi);

    case {'dp','dprime'}
        % D' = -(1/2pi) Im( n_x I_H(tau)(x) ).  Direct zk->0 limit of
        % helm2d_close 'dp': the Mlog piece vanishes (zk^2 prefactor),
        % leaving M = -Re(n_x .* HypC) / (2*pi).  HypC already has
        % source weights baked in via wlchs_target's thyp.
        if nargin < 8 || isempty(HypC)
            error('CHNK.HELSINGO.LAP2D_CLOSE: ''dp'' requires HypC');
        end
        if nargin < 9 || isempty(nzt)
            error('CHNK.HELSINGO.LAP2D_CLOSE: ''dp'' requires target normals nzt');
        end
        nzt = nzt(:);
        M = -real(nzt .* HypC) / (2*pi);

    otherwise
        error('CHNK.HELSINGO.LAP2D_CLOSE: unsupported type %s', type);
end
end
