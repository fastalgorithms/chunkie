function [mat] = kern(lam,mu,s,t,type)
%ELAST2D.KERN
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
% - [Dalt]_ij = -zeta*(r_ir_j r_n)/r^4 - 2*eta*r_n/r^4
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
mat = zeros(2*m,2*n);

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

end
