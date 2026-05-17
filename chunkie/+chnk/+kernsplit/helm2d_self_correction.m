function [delta, U_src_cell] = helm2d_self_correction(chnkr, type, zk, U_src_cell)
%CHNK.HELSINGO.HELM2D_SELF_CORRECTION  self-panel close-correction delta.
%
% Returns a sparse matrix `delta` such that the self-panel entries are the
% kernel-split correction to add to smooth Gauss-Legendre quadrature.  This
% ports testhelmos's SoperC_targ/ToperC_targ self-panel formulas to
% chunkie's G = (i/4) H_0 convention.
%
% Supported types: 's', 'sp', 'dp'.  'd' returns zero (smooth GL is exact
% on a smooth curve for D's diagonal limit).  Other types return zero too.

ngl = chnkr.k;
nch = chnkr.nch;
npt = chnkr.npt;

tlow = lower(type);
if ~ismember(tlow, {'s','single','sp','sprime','dp','dprime','t'})
    delta = sparse(npt, npt);
    return
end

[T, wgl] = lege.exps(ngl);
[WfrakL, accept] = chnk.kernsplit.wfrak_log(ngl, 0, 1);
if numel(accept) ~= ngl
    error('helm2d_self_correction: self log moments should accept all GL nodes');
end
LogC0 = -log(abs(T(:) - T(:).'));
LogC0(1:ngl+1:end) = 0;
LogC0 = LogC0 + WfrakL ./ wgl(:).';

[rends, ~] = chunkends(chnkr);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));

if nargin < 4 || isempty(U_src_cell)
    U_src_cell = cell(1, nch);
end

gamma = 0.5772156649015328606;
ii = [];
jj = [];
vv = [];

for ip = 1:nch
    idx = (ip-1)*ngl + (1:ngl);
    zsrc = chnkr.r(1,:,ip) + 1i*chnkr.r(2,:,ip);
    nz = chnkr.n(1,:,ip) + 1i*chnkr.n(2,:,ip);
    dz = chnkr.d(1,:,ip) + 1i*chnkr.d(2,:,ip);
    wzp = dz .* wgl(:).';
    awzp = chnkr.wts(:,ip).';
    logspd = log(abs(dz));

    ztgt = zsrc(:);
    nzt = nz(:);
    zdiff = ztgt - zsrc;

    switch tlow
        case {'s','single'}
            Mblk = -besselj(0, zk*abs(zdiff)) .* LogC0 / (2*pi);
            Mblk = Mblk .* awzp;

            tmp = 0.25i - (log(zk/2) + gamma) / (2*pi);
            dval = (tmp - (diag(LogC0).' + logspd) / (2*pi)) .* awzp;

        case {'sp','sprime'}
            % Helsing demo12.m / demo13b.m Kadjoperinit, scaled by 1/2 for
            % chunkie's G = (i/4) H_0 convention.
            % Off-diagonal: smooth coefficient of log r is
            %   (zk r) J_1(zk r) real(nzt/zdiff) awzp.
            % Use besselj(1,·) directly so the formula handles complex zk.
            tmp = zk * abs(zdiff);
            cosNT = real(nzt ./ zdiff);     % NaN on the diagonal i=j
            cosNT(1:ngl+1:end) = 0;          % clear diagonal: dval overwrites
            Mblk = (tmp .* besselj(1,tmp) .* cosNT) .* LogC0 .* awzp / (2*pi);
            % Diagonal: bare-kernel curvature limit S'(z,z) = -kappa/(4 pi).
            % awzp(i)*S'(z_i, z_i) = -imag(zpp/zp)(i) * w(i) / (4 pi).
            % chnkr.d2 -> z'' (canonical second derivative).
            zpp = chnkr.d2(1,:,ip) + 1i*chnkr.d2(2,:,ip);
            dval = -imag(zpp ./ dz) .* wgl(:).' / (4*pi);

        case {'dp','dprime','t'}
            tmp = zk * abs(zdiff);
            Ta1 = zk^2 * besselj(1,tmp) ./ tmp .* real(nzt * conj(nz));
            D1 = -real(nz ./ zdiff);
            D2 =  real(nzt ./ zdiff);
            Tb1 = tmp.^2 .* besselj(2,tmp) .* D1 .* D2;
            Mblk = -(Ta1 + Tb1) .* LogC0 / (2*pi);
            Mblk = Mblk .* awzp;

            if isempty(U_src_cell{ip})
                U_src_cell{ip} = chnk.kernsplit.wlchs_src_precomp( ...
                    endsa(ip), endsb(ip), zsrc);
            end
            U_src = U_src_cell{ip};
            [~, ~, HypC] = chnk.kernsplit.wlchs_target( ...
                endsa(ip), endsb(ip), ztgt, zsrc, nz, wzp, awzp, U_src, 1);
            Mblk = Mblk - real(nzt .* HypC) / (2*pi);

            tmpdiag = 0.5i - (log(zk/2) + gamma - 0.5) / pi;
            dval = zk^2 * 0.25 * ...
                (tmpdiag - (diag(LogC0).' + logspd) / pi) .* awzp;
            dval = dval - diag(real(nzt .* HypC) / (2*pi)).';
    end

    Mblk(1:ngl+1:end) = dval;
    [I, J] = ndgrid(idx, idx);
    ii = [ii; I(:)];
    jj = [jj; J(:)];
    vv = [vv; Mblk(:)];
end

delta = sparse(ii, jj, vv, npt, npt);
end
