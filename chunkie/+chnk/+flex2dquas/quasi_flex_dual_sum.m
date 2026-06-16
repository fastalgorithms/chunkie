function [val, grad, hess, third, fourth] = quasi_flex_dual_sum(rx,ry,zk,kappa,d)
%CHNK.FLEX2DQUAS.QUASI_FLEX_DUAL_SUM evaluate the quasiperiodic dual
% (Ewald) sum for the flexural wave Green's function using plane-wave
% expansions.
%
% Computes the far-field part of the quasiperiodic flexural Green's function
% as the difference of two quasiperiodic Helmholtz dual sums evaluated at
% wavenumbers zk and i*zk:
%
%   G_flex_far(x,y) = G_H_far(x,y; zk) - G_H_far(x,y; i*zk)
%
% This is used by CHNK.FLEX2DQUAS.LATTICECOEFS to determine the lattice
% sum coefficients.
%
% Syntax: [val,grad,hess,third,fourth] = ...
%             chnk.flex2dquas.quasi_flex_dual_sum(rx,ry,zk,kappa,d)
%
% Input:
%   rx    - x-components of evaluation points (targets relative to source)
%   ry    - y-components of evaluation points; must be nonzero
%   zk    - flexural wavenumber
%   kappa - (nkappa,1) quasiperiodic phase parameters
%   d     - period (scalar)
%
% Output:
%   val   - (nkappa, npt) Green's function values
%   grad  - (nkappa, npt, 2) gradient [d/drx, d/dry]
%   hess  - (nkappa, npt, 3) Hessian [d^2/drx^2, d^2/drxdry, d^2/dry^2]
%   third - (nkappa, npt, 4) third derivatives [xxx, xxy, xyy, yyy]
%   fourth- (nkappa, npt, 5) fourth derivatives [xxxx, ..., yyyy]
%
% see also CHNK.LAP2DQUAS.QUASI_DUAL_SUM, CHNK.FLEX2DQUAS.LATTICECOEFS

[valh0, gradh0, hessh0, thirdh0, fourthh0] = chnk.lap2dquas.quasi_dual_sum(rx,ry,zk,kappa,d);
[valk0, gradk0, hessk0, thirdk0, fourthk0] = chnk.lap2dquas.quasi_dual_sum(rx,ry,1i*zk,kappa,d);

over2k2 = 1; %/(2*zk.^2);
val = over2k2*(valh0-valk0);
grad = over2k2*(gradh0-gradk0);
hess = over2k2*(hessh0-hessk0);
third = over2k2*(thirdh0 - thirdk0);
fourth = over2k2*(fourthh0-fourthk0);

end