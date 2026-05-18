function ubdry = helmos_planewave_rhs(srcinfo, bctype, zk, theta)
%HELMOS_PLANEWAVE_RHS  Bruno-Lintner-compatible boundary data for an
% incident planewave u_inc = exp(i k d.r) with d = (cos theta, sin theta).
%
% For Dirichlet: u_total = 0 on Gamma -> u_scat = -u_inc -> ubdry = -u_inc.
% For Neumann:   du_total/dn = 0 -> du_scat/dn = -du_inc/dn
%                 -> ubdry = -i k (d.n) exp(i k d.r).
%
% Inputs: srcinfo (.r 2xN, .n 2xN), bctype 'd'/'n', zk, theta (rad).

x = srcinfo.r(1,:).';  y = srcinfo.r(2,:).';
d = [cos(theta); sin(theta)];

phase = exp(1i*abs(zk)*(x*d(1) + y*d(2)));
if bctype == 'd'
    ubdry = -phase;
elseif bctype == 'n'
    nx = srcinfo.n(1,:).';  ny = srcinfo.n(2,:).';
    ubdry = -1i*zk*phase.*(nx*d(1) + ny*d(2));
else
    error('helmos_planewave_rhs: bctype must be ''d'' or ''n''');
end
end
