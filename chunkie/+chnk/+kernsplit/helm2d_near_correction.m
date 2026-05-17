function delta = helm2d_near_correction(chnkr_src, chnkr_tgt, type, zk, dlim, coefs)
%CHNK.HELSINGO.HELM2D_NEAR_CORRECTION  off-panel close-correction delta.
%
% Returns a sparse ntgt-by-nsrc correction matrix to add to smooth
% Gauss-Legendre quadrature for target nodes on chnkr_tgt and source panels
% on chnkr_src. This ports testhelmos's osbuildmat_nearcorrection path to
% chunkie's G = (i/4)H_0 convention.
%
% Supported types: 's' / 'single', 'd' / 'double', 'sp' / 'sprime',
% 'dp' / 'dprime' / 't', 'c' / 'combined' (coefs(1)*D + coefs(2)*S),
% 'sc' / 'spcombined' (coefs(1)*S' + coefs(2)*S).

if nargin < 5 || isempty(dlim)
    dlim = 0.5;
end
if nargin < 6
    coefs = [];
end
if (strcmpi(type,'c') || strcmpi(type,'combined')) && (isempty(coefs) || numel(coefs) ~= 2)
    error("helm2d_near_correction: type 'c' requires coefs = [c1, c2]");
end
if (strcmpi(type,'sc') || strcmpi(type,'spcombined')) && (isempty(coefs) || numel(coefs) ~= 2)
    error("helm2d_near_correction: type 'sc' requires coefs = [c_sp, c_s]");
end

ngl = chnkr_src.k;
nch_src = chnkr_src.nch;
nsrc = chnkr_src.npt;
ntgt = chnkr_tgt.npt;

[~, wgl] = lege.exps(ngl);

ztgt_all = chnkr_tgt.r(1,:) + 1i*chnkr_tgt.r(2,:);
nzt_all = chnkr_tgt.n(1,:) + 1i*chnkr_tgt.n(2,:);

[rends, ~] = chunkends(chnkr_src);
endsa = squeeze(rends(1,1,:) + 1i*rends(2,1,:));
endsb = squeeze(rends(1,2,:) + 1i*rends(2,2,:));

ii = [];
jj = [];
vv = [];

for ip = 1:nch_src
    src_idx = (ip-1)*ngl + (1:ngl);
    zsrc = chnkr_src.r(1,:,ip) + 1i*chnkr_src.r(2,:,ip);
    nz = chnkr_src.n(1,:,ip) + 1i*chnkr_src.n(2,:,ip);
    dz = chnkr_src.d(1,:,ip) + 1i*chnkr_src.d(2,:,ip);
    wzp = dz .* wgl(:).';
    awzp = chnkr_src.wts(:,ip).';

    panlen = sum(awzp);
    zdiff = ztgt_all(:) - zsrc;
    d2 = real(zdiff).^2 + imag(zdiff).^2;
    near = min(d2, [], 2) < (dlim*panlen)^2;
    if ~any(near)
        continue
    end

    tgt_idx = find(near).';
    ztgt = ztgt_all(tgt_idx).';
    U = chnk.kernsplit.wlchs_src_precomp(endsa(ip), endsb(ip), zsrc);

    switch lower(type)
        case {'s','single'}
            LogC = chnk.kernsplit.wlchs_target( ...
                endsa(ip), endsb(ip), ztgt, zsrc, nz, wzp, awzp, U, 0);
            Mc = chnk.kernsplit.helm2d_close('s', zk, ztgt, zsrc, nz, wzp, LogC, []);
            Mblk = Mc .* awzp;

        case {'d','double'}
            [LogC, CauC] = chnk.kernsplit.wlchs_target( ...
                endsa(ip), endsb(ip), ztgt, zsrc, nz, wzp, awzp, U, 0);
            Mblk = chnk.kernsplit.helm2d_close('d', zk, ztgt, zsrc, nz, wzp, LogC, CauC);

        case {'sp','sprime'}
            % S' near correction: LogC + CauC (Helsing-Karlsson 2018 §4.4).
            [LogC, CauC] = chnk.kernsplit.wlchs_target( ...
                endsa(ip), endsb(ip), ztgt, zsrc, nz, wzp, awzp, U, 0);
            nzt = nzt_all(tgt_idx).';
            Mblk = chnk.kernsplit.helm2d_close('sp', zk, ztgt, zsrc, nz, ...
                                               wzp, LogC, CauC, [], nzt, awzp);

        case {'dp','dprime','t'}
            [LogC, ~, HypC] = chnk.kernsplit.wlchs_target( ...
                endsa(ip), endsb(ip), ztgt, zsrc, nz, wzp, awzp, U, 0);
            nzt = nzt_all(tgt_idx).';
            Mblk = chnk.kernsplit.helm2d_close('dp', zk, ztgt, zsrc, nz, ...
                                               wzp, LogC, [], HypC, nzt, awzp);

        case {'c','combined'}
            % combined kernel = coefs(1)*D + coefs(2)*S; near correction is
            % the linear superposition of D and S near corrections (kernel
            % split is linear in the kernel).
            [LogC, CauC] = chnk.kernsplit.wlchs_target( ...
                endsa(ip), endsb(ip), ztgt, zsrc, nz, wzp, awzp, U, 0);
            Mblk_d = chnk.kernsplit.helm2d_close('d', zk, ztgt, zsrc, nz, wzp, LogC, CauC);
            Mblk_s = chnk.kernsplit.helm2d_close('s', zk, ztgt, zsrc, nz, wzp, LogC, []) .* awzp;
            Mblk = coefs(1) * Mblk_d + coefs(2) * Mblk_s;

        case {'sc','spcombined'}
            % combined kernel = coefs(1)*S' + coefs(2)*S; same linearity.
            [LogC, CauC] = chnk.kernsplit.wlchs_target( ...
                endsa(ip), endsb(ip), ztgt, zsrc, nz, wzp, awzp, U, 0);
            nzt = nzt_all(tgt_idx).';
            Mblk_sp = chnk.kernsplit.helm2d_close('sp', zk, ztgt, zsrc, nz, ...
                                                  wzp, LogC, CauC, [], nzt, awzp);
            Mblk_s  = chnk.kernsplit.helm2d_close('s', zk, ztgt, zsrc, nz, wzp, LogC, []) .* awzp;
            Mblk = coefs(1) * Mblk_sp + coefs(2) * Mblk_s;

        otherwise
            error('helm2d_near_correction: unsupported type %s', type);
    end

    [I, J] = ndgrid(tgt_idx, src_idx);
    ii = [ii; I(:)];
    jj = [jj; J(:)];
    vv = [vv; Mblk(:)];
end

delta = sparse(ii, jj, vv, ntgt, nsrc);
end
