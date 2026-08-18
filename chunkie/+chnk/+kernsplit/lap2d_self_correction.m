function [delta, U_src_cell] = lap2d_self_correction(chnkr, type, U_src_cell)
%CHNK.HELSINGO.LAP2D_SELF_CORRECTION  self-panel close-correction for 2D
% Laplace layer potentials (Helsing-Ojala / wLCHS).
%
% Returns a sparse matrix `delta` whose diagonal entries are the FULL
% close-evaluated matrix entry, and off-diagonal entries are the kernel
% delta over smooth Gauss-Legendre quadrature.  The caller (chunkermat
% wlchs path) supplies M_smooth with diag zeroed; M_smooth + delta is
% then the full close-evaluated self block.
%
% Mirrors the zk -> 0 limit of helm2d_self_correction.

ngl = chnkr.k;
nch = chnkr.nch;
npt = chnkr.npt;
tlow = lower(type);
if ~ismember(tlow, {'s','single','d','double','sp','sprime','dp','dprime'})
    delta = sparse(npt, npt);
    return
end

if nargin < 3 || isempty(U_src_cell), U_src_cell = cell(1, nch); end

[T, wgl] = lege.exps(ngl);
[WfrakL, accept] = chnk.kernsplit.wfrak_log(ngl, 0, 1);
if numel(accept) ~= ngl
    error('lap2d_self_correction: self log moments should accept all GL nodes');
end
LogC0 = -log(abs(T(:) - T(:).'));
LogC0(1:ngl+1:end) = 0;
LogC0 = LogC0 + WfrakL ./ wgl(:).';

[rends, ~] = chunkends(chnkr);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));

ii = []; jj = []; vv = [];
for ip = 1:nch
    idx = (ip-1)*ngl + (1:ngl);
    zsrc = chnkr.r(1,:,ip) + 1i*chnkr.r(2,:,ip);
    nz   = chnkr.n(1,:,ip) + 1i*chnkr.n(2,:,ip);
    dz   = chnkr.d(1,:,ip) + 1i*chnkr.d(2,:,ip);
    awzp = chnkr.wts(:,ip).';
    wzp  = dz .* wgl(:).';
    logspd = log(abs(dz));
    a = endsa(ip); b = endsb(ip);

    ztgt = zsrc(:); nzt = nz(:);

    switch tlow
        case {'s','single'}
            % off-diag: M = -LogC0/(2*pi) * awzp  (matches zk->0 of helm)
            % diag:    M = -(log|dz| + diag(LogC0))/(2*pi) * awzp
            Mblk = -LogC0/(2*pi);
            Mblk = Mblk .* awzp;
            dval = -(logspd + diag(LogC0).')/(2*pi) .* awzp;
            Mblk(1:ngl+1:end) = dval;

        case {'d','double'}
            % Laplace D self block: needs CauC moment (wlchs_target self)
            if isempty(U_src_cell{ip})
                U_src_cell{ip} = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc);
            end
            U = U_src_cell{ip};
            [~, CauC] = chnk.kernsplit.wlchs_target(a, b, ztgt, zsrc, ...
                nz, wzp, awzp, U, 1);
            Mblk = chnk.kernsplit.lap2d_close('d', ztgt, zsrc, nz, wzp, ...
                LogC0, CauC);

        case {'sp','sprime'}
            if isempty(U_src_cell{ip})
                U_src_cell{ip} = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc);
            end
            U = U_src_cell{ip};
            [~, CauC] = chnk.kernsplit.wlchs_target(a, b, ztgt, zsrc, ...
                nz, wzp, awzp, U, 1);
            Mblk = chnk.kernsplit.lap2d_close('sp', ztgt, zsrc, nz, wzp, ...
                LogC0, CauC, [], nzt);

        case {'dp','dprime'}
            if isempty(U_src_cell{ip})
                U_src_cell{ip} = chnk.kernsplit.wlchs_src_precomp(a, b, zsrc);
            end
            U = U_src_cell{ip};
            [~, ~, HypC] = chnk.kernsplit.wlchs_target(a, b, ztgt, zsrc, ...
                nz, wzp, awzp, U, 1);
            Mblk = chnk.kernsplit.lap2d_close('dp', ztgt, zsrc, nz, wzp, ...
                LogC0, [], HypC, nzt);
    end

    [I, J] = ndgrid(idx, idx);
    ii = [ii; I(:)]; jj = [jj; J(:)]; vv = [vv; Mblk(:)];
end
delta = sparse(ii, jj, vv, npt, npt);
end
