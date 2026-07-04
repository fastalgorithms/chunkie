function submat = kern(prm, srcinfo, targinfo, type, varargin)
%CHNK.LNSF2D.KERN kernels for the 2D time-harmonic linearized
% Navier-Stokes-Fourier (thermoviscous acoustic) system, eliminated to the
% (v1,v2,T) form (time convention exp(-i omega t)):
%
%   mu Lap v + (lam+mu) grad div v + i om rho0 v - eta_v grad T = 0
%   kap Lap T + i om rho0 cv T - eta_T div v = 0
%
% The fundamental matrix (z = x-y, Gk = (i/4) H_0^(1)(k|z|)) is
%
%   Gvv = (1/mu) Gkv I + hess Phi,  Phi = (Gkv - ba Gka - bh Gkh)/(i om rho0)
%   GvT = -cvT grad(Gka - Gkh),     GTv = -cTv grad(Gka - Gkh)^T
%   GTT = (Ba Gka + Bh Gkh)/kap
%
% and the double layer is built with the Kupradze-type PSEUDOSTRESS
% conormal, theta = mu(lam+mu)/(lam+3mu), acting at the source point on the
% density index, which cancels the Cauchy-singular part of the kernel so
% that all blocks are weakly singular (log) on smooth boundaries:
%
%   Wvv_ij = -[mu n_k d_k Gvv_ij + th n_k d_j Gvv_ik + (lam+mu-th) n_j d_l Gvv_il]
%   WvT_i  = -kap n_k d_k GvT_i
%   WTv_j  = -[mu n_k d_k GTv_j + th n_k d_j GTv_k + (lam+mu-th) n_j d_l GTv_l]
%   WTT    = -kap n_k d_k GTT
%
% (all derivatives in z = x - y; n = n(y), the normal at the source point,
% pointing out of the obstacle). Stable near-diagonal evaluation uses
% CHNK.HELM2D.HELMDIFFGREEN for each of the three wavenumbers; the log
% parts recombine exactly because ba+bh = 1 and Ba+Bh = 1.
%
% Exterior Dirichlet (no-slip + isothermal) scattering: with the combined
% kernel C = D + eta S, the system is (1/2) I + chunkermat(chnkr, C), a
% second-kind system, and the representation is u = (D + eta S) psi.
%
% Syntax: submat = chnk.lnsf2d.kern(prm, srcinfo, targinfo, type, eta)
%
% prm - parameter struct from chnk.lnsf2d.params
% srcinfo, targinfo - ptinfo structs with fields r (2,n), and n (2,ns)
%                     required at sources for 'd' and 'c'
% type - 's' (single layer), 'd' (pseudostress double layer),
%        'c' (combined, D + eta*S, eta = varargin{1}, default 1)
%
% Output is (3 nt) x (3 ns), interleaved per node as (v1, v2, T).
%
% see also CHNK.LNSF2D.PARAMS, CHNK.HELM2D.HELMDIFFGREEN

src = srcinfo.r(1:2,:); targ = targinfo.r(1:2,:);
[~, ns] = size(src); [~, nt] = size(targ);

rx = targ(1,:).' - src(1,:);
ry = targ(2,:).' - src(2,:);
r2 = rx.^2 + ry.^2;
L0 = -log(r2)/(4*pi);          % -log(r)/(2 pi)

mu = prm.mu; lam = prm.lam; kap = prm.kap; th = prm.theta;
lmt = lam + mu - th;
c0 = prm.c0; ba = prm.ba; bh = prm.bh; Ba = prm.Ba; Bh = prm.Bh;

switch lower(type)

  case {'s', 'single'}

    [~, ga, ha] = chnk.helm2d.helmdiffgreen(prm.ka, src, targ);
    [~, gh, hh] = chnk.helm2d.helmdiffgreen(prm.kh, src, targ);
    [~, ~, hv] = chnk.helm2d.helmdiffgreen(prm.kv, src, targ);

    % NB: Gkv/mu and GTT are evaluated directly via H0 rather than through
    % the split L0 + helmdiffgreen value: beyond the boundary layer the
    % log parts cancel to exponentially small values, and the rounding
    % noise eps*|log r| would be amplified by 1/mu (resp. 1/kap) --
    % an O(eps*|log r|/mu) error that dominates manufactured-solution
    % evaluations at small viscosity.
    r = sqrt(r2);
    Gkvmu = (1i/4)*besselh(0, 1, prm.kv*r)/mu;
    Hc1 = c0*(hv(:,:,1) - ba*ha(:,:,1) - bh*hh(:,:,1));
    Hc2 = c0*(hv(:,:,2) - ba*ha(:,:,2) - bh*hh(:,:,2));
    Hc3 = c0*(hv(:,:,3) - ba*ha(:,:,3) - bh*hh(:,:,3));
    Gvv11 = Gkvmu + Hc1;
    Gvv12 = Hc2;
    Gvv22 = Gkvmu + Hc3;
    dgx = ga(:,:,1) - gh(:,:,1);
    dgy = ga(:,:,2) - gh(:,:,2);
    % beyond the boundary layer the two GD gradients share an O(1/r)
    % static part that cancels in the difference; the rounding left
    % over is amplified by cTv ~ 1/kap in the T rows. There compute
    % grad(G_a - G_h) directly from H1 instead (no static part). Near
    % the boundary layer the GD form is the accurate one.
    ifar = abs(prm.kh*r) > 4.0;
    if any(ifar(:))
        chfac = -(1i/4)*(prm.ka*besselh(1, 1, prm.ka*r(ifar)) ...
                       - prm.kh*besselh(1, 1, prm.kh*r(ifar)))./r(ifar);
        dgx(ifar) = chfac.*rx(ifar);
        dgy(ifar) = chfac.*ry(ifar);
        % the v rows of the split (helmdiffgreen) form inherit
        % ~1e-12..1e-10 relative noise from the GD evaluations; beyond
        % the boundary layer build Gvv directly from Hankel functions
        % (verified against an extended-precision reference to 3e-16)
        ra = r(ifar);
        H0a = besselh(0, 1, prm.ka*ra); H1a = besselh(1, 1, prm.ka*ra);
        H0h = besselh(0, 1, prm.kh*ra); H1h = besselh(1, 1, prm.kh*ra);
        H0v = besselh(0, 1, prm.kv*ra); H1v = besselh(1, 1, prm.kv*ra);
        G1a = -(1i/4)*prm.ka*H1a;
        G2a = -(1i/4)*prm.ka^2*(H0a - H1a./(prm.ka*ra));
        G1h = -(1i/4)*prm.kh*H1h;
        G2h = -(1i/4)*prm.kh^2*(H0h - H1h./(prm.kh*ra));
        G1v = -(1i/4)*prm.kv*H1v;
        G2v = -(1i/4)*prm.kv^2*(H0v - H1v./(prm.kv*ra));
        Phi1 = c0*(G1v - ba*G1a - bh*G1h);
        Phi2 = c0*(G2v - ba*G2a - bh*G2h);
        P1f = Phi1./ra; Q1f = Phi2 - P1f;
        Gv0mu = (1i/4)*H0v/mu;
        zxu = rx(ifar)./ra; zyu = ry(ifar)./ra;
        Gvv11(ifar) = Gv0mu + P1f + zxu.^2.*Q1f;
        Gvv22(ifar) = Gv0mu + P1f + zyu.^2.*Q1f;
        Gvv12(ifar) = zxu.*zyu.*Q1f;
    end
    GvT1 = -prm.cvT*dgx; GvT2 = -prm.cvT*dgy;
    GTv1 = -prm.cTv*dgx; GTv2 = -prm.cTv*dgy;
    GTT = (1i/4)*(Ba*besselh(0, 1, prm.ka*r) ...
                + Bh*besselh(0, 1, prm.kh*r))/kap;

    submat = zeros(3*nt, 3*ns);
    submat(1:3:end, 1:3:end) = Gvv11; submat(1:3:end, 2:3:end) = Gvv12;
    submat(2:3:end, 1:3:end) = Gvv12; submat(2:3:end, 2:3:end) = Gvv22;
    submat(1:3:end, 3:3:end) = GvT1;  submat(2:3:end, 3:3:end) = GvT2;
    submat(3:3:end, 1:3:end) = GTv1;  submat(3:3:end, 2:3:end) = GTv2;
    submat(3:3:end, 3:3:end) = GTT;

  case {'d', 'double', 'dps'}

    [~, ga, ha, ta] = chnk.helm2d.helmdiffgreen(prm.ka, src, targ);
    [~, gh, hh, th3] = chnk.helm2d.helmdiffgreen(prm.kh, src, targ);
    [~, ~, ~, tv] = chnk.helm2d.helmdiffgreen(prm.kv, src, targ);

    nx = repmat(srcinfo.n(1,:), nt, 1);
    ny = repmat(srcinfo.n(2,:), nt, 1);

    r = sqrt(r2);
    gL0x = -rx./(2*pi*r2); gL0y = -ry./(2*pi*r2);
    % grad G_kv directly from H1 (no log split): the split form
    % gL0 + grad GD_kv cancels two O(1/r)-scale parts to exponentially
    % small beyond the boundary layer, and the rounding left over,
    % divided by mu, is a delta^{-2}-scaled noise floor in the v rows
    gkfac = -(1i/4)*prm.kv*besselh(1, 1, prm.kv*r)./r;
    gkv = cat(3, gkfac.*rx, gkfac.*ry);
    % symmetric 3rd-derivative combo, storage (xxx,xxy,xyy,yyy)
    T3c = c0*(tv - ba*ta - bh*th3);
    % dGvv(i,j,k) = d_k Gvv_ij = (1/mu) gkv_k delta_ij + T3c(i,j,k)
    d3idx = @(i,j,k) (i==2) + (j==2) + (k==2) + 1;
    dGvv = @(i,j,k) T3c(:,:,d3idx(i,j,k)) + (i==j)*gkv(:,:,k)/mu;
    % trace: sum_l d_l Gvv_il = c0(ba ka^2 Ga' + bh kh^2 Gh')(z_i/r)
    % exactly, by Lap G_k = -k^2 G_k and mu c0 kv^2 = 1. The assembled
    % form subtracts two 1/mu-scale quantities in floating point.
    gplfac = c0*(ba*prm.ka^2*(-(1i/4)*prm.ka*besselh(1, 1, prm.ka*r)) ...
               + bh*prm.kh^2*(-(1i/4)*prm.kh*besselh(1, 1, prm.kh*r)))./r;

    Wvv = cell(2,2);
    for i = 1:2
      for j = 1:2
        t1 = mu*(nx.*dGvv(i,j,1) + ny.*dGvv(i,j,2));
        t2 = th*(nx.*dGvv(i,1,j) + ny.*dGvv(i,2,j));
        nj = nx*(j==1) + ny*(j==2);
        zi = rx*(i==1) + ry*(i==2);
        t3 = lmt*nj.*(gplfac.*zi);
        Wvv{i,j} = -(t1 + t2 + t3);
      end
    end

    % coupling hessian difference, storage (xx,xy,yy)
    dH = ha - hh;
    hidx = @(i,k) (i==2) + (k==2) + 1;
    HD = @(i,k) dH(:,:,hidx(i,k));
    WvT = cell(2,1); WTv = cell(2,1);
    trHD = HD(1,1) + HD(2,2);
    for i = 1:2
      WvT{i} = kap*prm.cvT*(nx.*HD(i,1) + ny.*HD(i,2));
    end
    for j = 1:2
      nj = nx*(j==1) + ny*(j==2);
      WTv{j} = prm.cTv*((mu+th)*(HD(j,1).*nx + HD(j,2).*ny) + lmt*nj.*trHD);
    end
    gTTx = gL0x + Ba*ga(:,:,1) + Bh*gh(:,:,1);
    gTTy = gL0y + Ba*ga(:,:,2) + Bh*gh(:,:,2);
    WTT = -(nx.*gTTx + ny.*gTTy);

    submat = zeros(3*nt, 3*ns);
    submat(1:3:end, 1:3:end) = Wvv{1,1}; submat(1:3:end, 2:3:end) = Wvv{1,2};
    submat(2:3:end, 1:3:end) = Wvv{2,1}; submat(2:3:end, 2:3:end) = Wvv{2,2};
    submat(1:3:end, 3:3:end) = WvT{1};   submat(2:3:end, 3:3:end) = WvT{2};
    submat(3:3:end, 1:3:end) = WTv{1};   submat(3:3:end, 2:3:end) = WTv{2};
    submat(3:3:end, 3:3:end) = WTT;

  case {'c', 'comb', 'combined'}

    % optional trailing arguments: eta (combined-field coupling, default 1)
    % and tscale (thermal-density conjugation, default 1). With tscale
    % the kernel is diag(1,1,1/tscale) * K * diag(1,1,tscale): the
    % (T,v) blocks scale like eta_T/kappa ~ delta^{-2} and the (v,T)
    % blocks like delta, so tscale ~ |kh| (dividing the T-rows by |kh|
    % and multiplying the T-columns) balances the system uniformly in
    % the boundary-layer thickness (measured: cond ~ 9, independent of
    % delta, versus cond ~ delta^{-2} unconjugated) (see the block-imbalance discussion in the
    % accompanying note). The conjugation commutes with the identity, so
    % jump relations and the rcip convention are unaffected; the physical
    % density is recovered as psi_T = tscale * psi_T'.
    eta = 1.0;
    if nargin > 4 && ~isempty(varargin{1}), eta = varargin{1}; end
    submat = chnk.lnsf2d.kern(prm, srcinfo, targinfo, 'd') + ...
             eta*chnk.lnsf2d.kern(prm, srcinfo, targinfo, 's');
    tscale = 1.0;
    if nargin > 5 && ~isempty(varargin{2}), tscale = varargin{2}; end
    if tscale ~= 1.0
        submat(3:3:end, :) = submat(3:3:end, :)/tscale;
        submat(:, 3:3:end) = submat(:, 3:3:end)*tscale;
    end

  otherwise
    error('CHNK.LNSF2D.KERN: unknown type %s', type);
end

end
