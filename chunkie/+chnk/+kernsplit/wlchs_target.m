function [wcmpL,wcmpC,wcmpH,wcmpS] = wlchs_target(a,b,ztg,zsc,nz,wzp,awzp,U,ifself)
%CHNK.HELSINGO.WLCHS_TARGET  panel-quadrature corrections (log, Cauchy,
% hypersingular, supersingular) for one panel and a vector of target
% points, in the kernel-split style of Helsing & Ojala.
%
% These corrections are MEANT TO BE ADDED to a smooth-Gauss-Legendre
% panel evaluation of the kernel.  For a target far from the panel the
% corrections vanish; for a target close to (or on) the panel they
% recover the full singular contribution.
%
% Inputs:
%   a, b   - complex endpoints of the panel (panel parameterized as a
%            complexified contour from a to b).
%   ztg    - nt x 1 complex target points.
%   zsc    - 1 x Ng complex source nodes (Gauss-Legendre on the panel).
%   nz     - 1 x Ng complex unit normals at source nodes.
%   wzp    - 1 x Ng complex contour weights z'(t) * w_GL.
%   awzp   - 1 x Ng real arclength weights |z'(t)| * w_GL.
%   U      - Ng x Ng inverse Vandermonde from wlchs_src_precomp.
%   ifself - 1 if the targets coincide with the source panel nodes
%            (self-interaction), 0 otherwise.
%
% Outputs (each of size nt x Ng, target by source):
%   wcmpL  - log-singular correction (always computed)
%   wcmpC  - Cauchy singular correction (computed if nargout >= 2)
%   wcmpH  - hypersingular correction  (computed if nargout >= 3)
%   wcmpS  - supersingular correction  (computed if nargout >= 4)
%
% O(p^2) per target.  Vectorized over targets.
%
% Ported from testhelmos.m / Bruno-Lintner-Helsing close-evaluation
% machinery.  See:
%   Helsing & Ojala, Corner singularities for elliptic problems...
%   J. Comput. Phys. 227 (2008) 8820-8840.

zsc  = zsc(:).';
nz   = nz(:).';
wzp  = wzp(:).';
awzp = awzp(:).';
ztg  = ztg(:);

Ng = numel(zsc);
nt = numel(ztg);
cc = (b-a)/2;
zt = (ztg - (b+a)/2)/cc;

% f(:,k+1) = \int_{-1}^{1} P_k(z)/(z - zt) dz, recurrence in zt
f = zeros(nt, Ng+1);
if ifself == 1
    upp = log(1 - zt);
    loo = log(-1 - zt);
    f(:,1) = upp - loo;
    ind = imag(f(:,1)) < -1.1;  f(ind,1) = f(ind,1) + 1i*pi;
    ind = imag(f(:,1)) >  1.1;  f(ind,1) = f(ind,1) - 1i*pi;
else
    % decide per target which side of the (complexified) panel it is on
    ifleft = false(nt,1);
    zdiff  = ztg - zsc;            % nt x Ng
    xd = real(zdiff); yd = imag(zdiff);
    d2 = xd.*xd + yd.*yd;
    [~, nearidx] = min(d2, [], 2);

    zsctr   = (zsc - (b+a)/2)/cc;
    imzsctr = imag(zsctr);
    [~, mxi] = max(abs(imzsctr));
    isconvex = sign(imzsctr(mxi)) == -1;

    nz_re = real(nz); nz_im = imag(nz);
    for k = 1:nt
        rezt = real(zt(k));
        imzt = imag(zt(k));
        nk = nearidx(k);
        costheta = xd(k,nk)*nz_re(nk) + yd(k,nk)*nz_im(nk);
        if costheta < 0
            ifleft(k) = true;
        end
        if abs(rezt) < 1
            if isconvex
                if costheta < 0 && imzt < 0 && imzt > 1.1*min(imzsctr)
                    if imzt < chnk.kernsplit.mydepth(zsctr,rezt,Ng)
                        ifleft(k) = false;
                    end
                end
            else
                if costheta > 0 && imzt > 0 && imzt < 1.1*max(imzsctr)
                    if imzt > chnk.kernsplit.mydepth(zsctr,rezt,Ng)
                        ifleft(k) = true;
                    end
                end
            end
        end
    end

    gam    = -1i*ones(nt,1);
    loggam = -0.5i*pi*ones(nt,1);
    gam(ifleft)    = 1i;
    loggam(ifleft) =  0.5i*pi;
    f(:,1) = loggam + log((1-zt)./(gam.*(-1-zt)));
end

f(:,2) = 2 + zt.*f(:,1);
for k = 1:Ng-1
    f(:,k+2) = ((2*k+1)*zt.*f(:,k+1) - k*f(:,k))/(k+1);
end

% q(:,k+1) = \int_{-1}^{1} P_k(z)*log(z - zt) dz
q = zeros(nt, Ng);
P111   = log(abs((1-zt).*(1+zt)));
q(:,1) = P111 - zt.*f(:,1) - 2;
q(:,2:Ng) = (f(:,1:Ng-1) - f(:,3:Ng+1))./(3:2:2*Ng-1);

% (nearly) log corrections
tlog = log(abs((zsc - ztg)/cc));   % nt x Ng
tlog(isinf(tlog)) = 0;
wcmpL = imag((q*U)*cc.*conj(nz))./awzp - tlog;

if nargout >= 2
    tcau = wzp./(zsc - ztg)/1i;
    tcau(isinf(tcau)) = 0;
    wcmpC = (f(:,1:Ng)*U)/1i - tcau;
end

if nargout >= 3
    fp = zeros(nt, Ng);
    fp(:,2) = f(:,1);
    for k = 1:Ng-2
        fp(:,k+2) = (2*k+1)*f(:,k+1) + fp(:,k);
    end
    c2 = -(1./(1-zt) + (-1).^(0:Ng-1)./(1+zt));
    fp = fp + c2;
    thyp = wzp./(zsc - ztg).^2/1i;
    thyp(isinf(thyp)) = 0;
    wcmpH = (fp*U)/(1i*cc) - thyp;
end

if nargout >= 4
    fpp = zeros(nt, Ng);
    for k = 1:Ng-2
        fpp(:,k+2) = (2*k+1)*fp(:,k+1) + fpp(:,k);
    end
    c3 = -0.5*(1./(1-zt).^2 + (-1).^(1:Ng)./(1+zt).^2) ...
         -0.25*(0:Ng-1).*(1:Ng).*(1./(1-zt) + (-1).^(1:Ng)./(1+zt));
    fpp = 0.5*fpp + c3;
    tsup = wzp./(zsc - ztg).^3/1i;
    tsup(isinf(tsup)) = 0;
    wcmpS = (fpp*U)/(1i*cc^2) - tsup;
end
end
