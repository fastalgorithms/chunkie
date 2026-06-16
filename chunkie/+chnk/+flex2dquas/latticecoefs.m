function S = latticecoefs(n,zk,d,kappa,alpha,a,M,l)
%CHNK.FLEX2DQUAS.LATTICECOEFS precompute lattice sum coefficients for the
% quasiperiodic flexural wave Green's function.
%
% The flexural Green's function is the difference G_H(zk) - G_H(i*zk) of
% two quasiperiodic Helmholtz Green's functions, so this routine calls
% CHNK.HELM2DQUAS.LATTICECOEFS twice and concatenates the results.
%
% Syntax: S = chnk.flex2dquas.latticecoefs(n,zk,d,kappa,alpha,a,M,l)
%
% Input:
%   n     - (N+1,1) vector of orders, typically 0:N
%   zk    - complex number, flexural wavenumber
%   d     - period (scalar)
%   kappa - (nkappa,1) array of quasiperiodic phase parameters
%   alpha - (nkappa,1) array exp(1i*kappa*d), Bloch factors
%   a     - upper limit of trapezoid quadrature interval (default 15)
%   M     - number of trapezoid quadrature nodes (default 1e4)
%   l     - number of explicit periodic copies on each side to exclude
%               (default 2)
%
% Output:
%   S - (nkappa, N+1, 2) lattice sum coefficients; S(:,:,1) are the
%       Helmholtz coefficients at wavenumber zk, S(:,:,2) at wavenumber i*zk
%
% see also CHNK.HELM2DQUAS.LATTICECOEFS, CHNK.FLEX2DQUAS.GREEN

sn1 = chnk.helm2dquas.latticecoefs(n,zk,d,kappa,alpha,a,M,l);
sn2 = chnk.helm2dquas.latticecoefs(n,1i*zk,d,kappa,alpha,a,M,l);
S = cat(3,sn1,sn2);

end