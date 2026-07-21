function [K] = kern(lam,mu,s,t,type,ifimage)
%CHNK.ELAST2D.KERN
% elasticity kernels in two dimensions
%
% notation
%
% - delta_ij = kronecker delta
% - r_i(m,l) = t.r(i,m) - s.r(i,l)
% - r^2 = r_ir_i
% - r = sqrt(r^2)
% - n_i = t.n(i,m) (for Strac)
% - n_i = s.n(i,l) (for D, Dalt)
% - r_n = r_i(m,l) t.n(i,m) (for Strac)
% - r_n = r_i(m,l) s.n(i,l) (for D, Dalt)
% - beta = (lam+3*mu)/(4*pi*mu*(lam+2*mu))
% - gamma = -(lam+mu)/(4*pi*mu*(lam+2*mu))
% - zeta = (lam+mu)/(pi*(lam+2*mu))
% - eta = mu/(2*pi*(lam+2*mu))
%
% kernels
%
% - [S]_ij = (beta*log(r) + gamma/2) delta_ij + gamma*r_ir_j/r^2
% - [T]_ijk = eta(delta_jl r_i/r^2 + delta_ik r_j/r^2 ...
%                - delta_ij r_k/r^2) + zeta r_ir_jr_k/r^4
% - [Strac]_ij = T_{ikj} n_k
% - [D]_ij = -T_{jki}n_k
% - [Dalt]_ij = -zeta*(r_ir_j r_n)/r^4 - 2*eta*r_n/r^2*delta_ij
%
% inputs
% - lam, mu: real floats, Lame parameters
% - s: ptinfo struct with source information
% - t: ptinfo struct with target information
% - type: string
%   * type == 'S', single layer kernel
%   * type == 'Strac', traction of single layer kernel
%   * type == 'Sgrad', gradient of single layer kernel
%                kernel is 4 x 2, where entries are organized as
%                          Sgrad = [grad S(1,1) grad S(1,2)
%                                   grad S(2,1) grad S(2,2)]
%   * type == 'D', double layer kernel
%   * type == 'Dalt', alternative (smooth) double layer kernel
%   * type == 'Daltgrad', gradient of alternative (smooth) double layer
%                kernel is 4 x 2, where entries are organized as
%                          Daltgrad = [grad Dalt(1,1) grad Dalt(1,2)
%                                      grad Dalt(2,1) grad Dalt(2,2)]
%   * type == 'Dalttrac', traction of alternative (smooth) double layer
%
% outpts
% - mat: (2 nt) x (2 ns) matrix
%

if nargin < 5
    type = 's';
end

if nargin < 6
    ifimage = false;
end

[~, ns] = size(s.r(:,:));
[~, nt] = size(t.r(:,:));

rx = t.r(1,:).' - s.r(1,:);
ry = t.r(2,:).' - s.r(2,:);

hp = zeros(1,ns);

if (ifimage)
    if (isfield(s,'data'))
        hp = s.data(1,:);
        if (size(s.data(:,:),1)>1)
            nx = s.data(2,:);
            ny = s.data(3,:);
        else
            nx = s.n(1,:);
            ny = s.n(2,:);
        end
    else
        error('missing data field in elastic strings');
    end
end

rimx = s.r(1,:) + hp.*nx;
rimy = s.r(2,:) + hp.*ny;

zx = t.r(1,:).' - rimx;
zy = t.r(2,:).' - rimy;

r2 = rx.^2 + ry.^2;
z2 = zx.^2 + zy.^2;

rfac = 1/2/pi/mu;


alph = (lam + mu)/(lam + 2*mu);

switch lower(type)
    case {'s'}
        nxs = nx;
        nys = ny;

        taux =  nys;
        tauy = -nxs;

        rt = rx.*taux + ry.*tauy;
        zt = zx.*taux + zy.*tauy;
        rn = rx.*nxs + ry.*nys;
        zn = zx.*nxs + zy.*nys;

        logr = -1/alph*log(r2)/2;
        logz = -1/alph*log(z2)/2;
        attr = atan2(rt, -rn)*(1-alph)/alph;
        attz = atan2(zt, -zn)*(1-alph)/alph;

        Kxx = rfac*(logr + rx.*rx./r2);
        Kyy = rfac*(logr + ry.*ry./r2);
        Kxy = rfac*(rx.*ry./r2 + attr);
        Kyx = rfac*(rx.*ry./r2 - attr);

        Wxx = rfac*(logz + zx.*zx./z2);
        Wyy = rfac*(logz + zy.*zy./z2);
        Wxy = rfac*(zx.*zy./z2 + attz);
        Wyx = rfac*(zx.*zy./z2 - attz);

    case {'strac'}
        rfuse = -4*rfac;
        nxt = t.n(1,:).';
        nyt = t.n(2,:).';
        rnt = rx.*nxt + ry.*nyt;
        znt = zx.*nxt + zy.*nyt;

        r4inv = 1./(r2.*r2);
        z4inv = 1./(z2.*z2);

        Kxx = rfuse.*(rx.*rx.*rnt).*r4inv;
        Kxy = rfuse.*(rx.*ry.*rnt).*r4inv;
        Kyx = rfuse.*(ry.*rx.*rnt).*r4inv;
        Kyy = rfuse.*(ry.*ry.*rnt).*r4inv;

        Wxx = rfuse.*(zx.*zx.*znt).*z4inv;
        Wxy = rfuse.*(zx.*zy.*znt).*z4inv;
        Wyx = rfuse.*(zy.*zx.*znt).*z4inv;
        Wyy = rfuse.*(zy.*zy.*znt).*z4inv;

    otherwise
        error('Unknown elasticity kernel type ''%s''.', type);
end


K = zeros(2*nt, 2*ns);
K(1:2:end, 1:2:end) = Kxx;
K(2:2:end, 1:2:end) = Kyx;
K(1:2:end, 2:2:end) = Kxy;
K(2:2:end, 2:2:end) = Kyy;

if (ifimage)
    K(1:2:end,1:2:end) = K(1:2:end,1:2:end) - Wxx;
    K(2:2:end,1:2:end) = K(2:2:end,1:2:end) - Wyx;
    K(1:2:end,2:2:end) = K(1:2:end,2:2:end) - Wxy;
    K(2:2:end,2:2:end) = K(2:2:end,2:2:end) - Wyy;
end

end
