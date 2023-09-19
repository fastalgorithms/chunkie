function [mat] = kern(lam,mu,s,t,type)
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

beta = (lam+3*mu)/(4*pi*mu*(lam+2*mu));
gamma = -(lam+mu)/(4*pi*mu*(lam+2*mu));
eta = mu/(2*pi*(lam+2*mu));
zeta = (lam+mu)/(pi*(lam+2*mu));

n = size(s.r,2);
m = size(t.r,2);

x = (t.r(1,:)).' - s.r(1,:);
y = (t.r(2,:)).' - s.r(2,:);

r2 = x.^2 + y.^2;
r4 = r2.^2;

if (strcmpi(type,'s'))
    logr = beta*log(r2)/2;
    rirjor2 = (kron(gamma*x./r2,[1 0; 1 0]) + ...
        kron(gamma*y./r2,[0 1; 0 1])).* ...
        (kron(x,[1 1; 0 0]) + kron(y,[0 0; 1 1]));
    mat = kron(logr+gamma/2,eye(2)) + rirjor2;
end
if (strcmpi(type,'strac'))
    dirx = t.n(1,:); dirx = dirx(:);
    diry = t.n(2,:); diry = diry(:);
    rdotv = dirx(:).*x + diry(:).*y;
    term2 = zeta*(kron(x.*rdotv./r4,[1 0; 1 0]) + ...
        kron(y.*rdotv./r4,[0 1; 0 1])).* ...
        (kron(x,[1 1; 0 0]) + kron(y,[0 0; 1 1]));
    term1 = eta*(kron(x.*diry./r2,[0 1; -1 0]) + ...
        kron(y.*dirx./r2,[0 -1; 1 0]) + ...
        kron(rdotv./r2,[1 0; 0 1]));
    mat = term1+term2;
end
if (strcmpi(type,'sgrad'))
    mat = zeros(4*m,2*n);
    tmp = beta*x./r2;
    mat(1:4:end,1:2:end) = tmp;
    mat(3:4:end,2:2:end) = tmp;
    tmp = beta*y./r2;
    mat(2:4:end,1:2:end) = tmp;
    mat(4:4:end,2:2:end) = tmp;
    
    % x der of x^2/ r^2
    tmp = gamma*(2*r2.*x - x.^2.*(2*x))./r4;
    mat(1:4:end,1:2:end) = mat(1:4:end,1:2:end) + tmp;
    % y der of x^2/ r^2
    tmp = gamma*(-x.^2.*(2*y))./r4;
    mat(2:4:end,1:2:end) = mat(2:4:end,1:2:end) + tmp;
    % x der of xy/ r^2
    tmp = gamma*(r2.*y-x.*y.*(2*x))./r4;
    mat(3:4:end,1:2:end) = mat(3:4:end,1:2:end) + tmp;
    mat(1:4:end,2:2:end) = mat(1:4:end,2:2:end) + tmp;
    % y der of xy/ r^2
    tmp = gamma*(r2.*x-x.*y.*(2*y))./r4;
    mat(4:4:end,1:2:end) = mat(4:4:end,1:2:end) + tmp;
    mat(2:4:end,2:2:end) = mat(2:4:end,2:2:end) + tmp;
    % x der of y^2/ r^2
    tmp = gamma*(- y.^2.*(2*x))./r4;
    mat(3:4:end,2:2:end) = mat(3:4:end,2:2:end) + tmp;
    % y der of y^2/ r^2
    tmp = gamma*(2*r2.*y-y.^2.*(2*y))./r4;
    mat(4:4:end,2:2:end) = mat(4:4:end,2:2:end) + tmp;
    
end
if (strcmpi(type,'d'))
    dirx = s.n(1,:); dirx = dirx(:).';
    diry = s.n(2,:); diry = diry(:).';
    rdotv = x.*dirx + y.*diry;
    term2 = zeta*(kron(x.*rdotv./r4,[1 0; 1 0]) + ...
        kron(y.*rdotv./r4,[0 1; 0 1])).* ...
        (kron(x,[1 1; 0 0]) + kron(y,[0 0; 1 1]));
    term1 = eta*(kron(x.*diry./r2,[0 -1; 1 0]) + ...
        kron(y.*dirx./r2,[0 1; -1 0]) + ...
        kron(rdotv./r2,[1 0; 0 1]));
    mat = -(term1+term2);
end

if (strcmpi(type,'dalt'))
    dirx = s.n(1,:); dirx = dirx(:).';
    diry = s.n(2,:); diry = diry(:).';
    rdotv = x.*dirx + y.*diry;
    term2 = zeta*(kron(x.*rdotv./r4,[1 0; 1 0]) + ...
        kron(y.*rdotv./r4,[0 1; 0 1])).* ...
        (kron(x,[1 1; 0 0]) + kron(y,[0 0; 1 1]));
    term1 = eta*(kron(rdotv./r2,[2 0; 0 2]));
    mat = -(term1+term2);
end

if (strcmpi(type,'daltgrad'))
    dirx = s.n(1,:); dirx = dirx(:).';
    diry = s.n(2,:); diry = diry(:).';
    rdotv = x.*dirx + y.*diry;
    r6 = r4.*r2;
    
    mat = zeros(4*m,2*n);
    
    % i = 1, j = 1, l = 1
    aij_xl = -zeta*(-4*x.^3.*rdotv./r6+(2*x.*rdotv+x.^2.*dirx)./r4) ...
        -2*eta*(dirx./r2 - 2*rdotv.*x./r4);
    mat(1:4:end,1:2:end) = aij_xl;

    % i = 1, j = 1, l = 2
    aij_xl = -zeta*(-4*x.^2.*y.*rdotv./r6+(x.^2.*diry)./r4) ...
        -2*eta*(diry./r2 - 2*rdotv.*y./r4);
    mat(2:4:end,1:2:end) = aij_xl;

    % i = 2, j = 1, l = 1
    aij_xl = -zeta*(-4*x.^2.*y.*rdotv./r6+(y.*rdotv+x.*y.*dirx)./r4);
    mat(3:4:end,1:2:end) = aij_xl;

    % i = 2, j = 1, l = 2
    aij_xl = -zeta*(-4*x.*y.^2.*rdotv./r6+(x.*rdotv+x.*y.*diry)./r4);
    mat(4:4:end,1:2:end) = aij_xl;

    % i = 1, j = 2, l = 1
    aij_xl = -zeta*(-4*x.^2.*y.*rdotv./r6+(y.*rdotv+x.*y.*dirx)./r4);
    mat(1:4:end,2:2:end) = aij_xl;

    % i = 1, j = 2, l = 2
    aij_xl = -zeta*(-4*x.*y.^2.*rdotv./r6+(x.*rdotv + x.*y.*diry)./r4);
    mat(2:4:end,2:2:end) = aij_xl;

    % i = 2, j = 2, l = 1
    aij_xl = -zeta*(-4*y.^2.*x.*rdotv./r6+(y.^2.*dirx)./r4) ...
        -2*eta*(dirx./r2 - 2*rdotv.*x./r4);
    mat(3:4:end,2:2:end) = aij_xl;

    % i = 2, j = 2, l = 2
    aij_xl = -zeta*(-4*y.^3.*rdotv./r6+(2*y.*rdotv+y.^2.*diry)./r4) ...
        -2*eta*(diry./r2 - 2*rdotv.*y./r4);
    mat(4:4:end,2:2:end) = aij_xl;

end

if (strcmpi(type,'dalttrac'))
    dirx = s.n(1,:); dirx = dirx(:).';
    diry = s.n(2,:); diry = diry(:).';
    rdotv = x.*dirx + y.*diry;
    r6 = r4.*r2;
    
    matg = zeros(4*m,2*n);
    
    % i = 1, j = 1, l = 1
    aij_xl = -zeta*(-4*x.^3.*rdotv./r6+(2*x.*rdotv+x.^2.*dirx)./r4) ...
        -2*eta*(dirx./r2 - 2*rdotv.*x./r4);
    matg(1:4:end,1:2:end) = aij_xl;

    % i = 1, j = 1, l = 2
    aij_xl = -zeta*(-4*x.^2.*y.*rdotv./r6+(x.^2.*diry)./r4) ...
        -2*eta*(diry./r2 - 2*rdotv.*y./r4);
    matg(2:4:end,1:2:end) = aij_xl;

    % i = 2, j = 1, l = 1
    aij_xl = -zeta*(-4*x.^2.*y.*rdotv./r6+(y.*rdotv+x.*y.*dirx)./r4);
    matg(3:4:end,1:2:end) = aij_xl;

    % i = 2, j = 1, l = 2
    aij_xl = -zeta*(-4*x.*y.^2.*rdotv./r6+(x.*rdotv+x.*y.*diry)./r4);
    matg(4:4:end,1:2:end) = aij_xl;

    % i = 1, j = 2, l = 1
    aij_xl = -zeta*(-4*x.^2.*y.*rdotv./r6+(y.*rdotv+x.*y.*dirx)./r4);
    matg(1:4:end,2:2:end) = aij_xl;

    % i = 1, j = 2, l = 2
    aij_xl = -zeta*(-4*x.*y.^2.*rdotv./r6+(x.*rdotv + x.*y.*diry)./r4);
    matg(2:4:end,2:2:end) = aij_xl;

    % i = 2, j = 2, l = 1
    aij_xl = -zeta*(-4*y.^2.*x.*rdotv./r6+(y.^2.*dirx)./r4) ...
        -2*eta*(dirx./r2 - 2*rdotv.*x./r4);
    matg(3:4:end,2:2:end) = aij_xl;

    % i = 2, j = 2, l = 2
    aij_xl = -zeta*(-4*y.^3.*rdotv./r6+(2*y.*rdotv+y.^2.*diry)./r4) ...
        -2*eta*(diry./r2 - 2*rdotv.*y./r4);
    matg(4:4:end,2:2:end) = aij_xl;

    mat = zeros(2*m,2*n);
    n1 = t.n(1,:); n1 = n1(:); n2 = t.n(2,:); n2 = n2(:);
    tmp = matg(1:4:end,:)+matg(4:4:end,:);
    mat(1:2:end,:) = lam*n1.*tmp;
    mat(2:2:end,:) = lam*n2.*tmp;
    tmp = mu*(matg(2:4:end,:) + matg(3:4:end,:));
    mat(1:2:end,:) = mat(1:2:end,:) + tmp.*n2 + 2*mu*matg(1:4:end,:).*n1;
    mat(2:2:end,:) = mat(2:2:end,:) + tmp.*n1 + 2*mu*matg(4:4:end,:).*n2;

end

end
