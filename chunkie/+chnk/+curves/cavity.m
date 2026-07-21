function [r,d,d2] = cavity(t,varargin)
%CHNK.CURVES.CAVITY return position, first and second derivatives of a
% cavity-shaped domain, given by a radially-deformed ellipse in polar
% coordinates
%
% Syntax: [r,d,d2] = chnk.curves.cavity(t,j)
%
% Input:
%   t - array of points (in [0,2*pi])
%
% Optional input:
%   j - real number, controls eccentricity and "depth" of the cavity.
%       should satisfy 0 <= j <~ 11.12 to avoid self-intersection
%       (default 10)
%
% Output:
%   r - 2 x numel(t) array of positions, r(:,i) = [x(t(i)); y(t(i))]
%   d - 2 x numel(t) array of t derivative of r
%   d2 - 2 x numel(t) array of second t derivative of r
%
% Examples:
%   [r,d,d2] = chnk.curves.cavity(t); % get default settings
%   [r,d,d2] = chnk.curves.cavity(t,j); % change j

j = 10;
if nargin > 1 && ~isempty(varargin{1})
    j = varargin{1};
end

if j > 11
    warning('chnk.curves.cavity:selfintersect', ...
        'cavity curve may self intersect for j > 11');
end

rx = (1-(j-1)/11)+0.03*(j-1);
ry = (1-(j-1)/11)+0.8*(j-1);

ct = cos(t);
st = sin(t);

xt = rx*ct + 20;
dxdt = -rx*st;
d2xdt2 = -rx*ct;
yt = ry*st;
dydt = ry*ct;
d2ydt2 = -ry*st;

rmax = 400*rx*rx/(ry*ry-rx*rx) + 400 + ry^2;

a = (1+ 0.7*(j-1))/2;
C = 1/rmax^a;
rho_t = C*(xt.^2 + yt.^2).^(a);
drho_t_dt = 2*C*a*(xt.^2 + yt.^2).^(a-1).*(xt.*dxdt + yt.*dydt);
d2rho_t_dt2 = C*((a-1)*a*(xt.^2+yt.^2).^(a-2).*(2*xt.*dxdt + 2*yt.*dydt).^2 + ...
     a*(xt.*xt + yt.*yt).^(a-1).*(2*xt.*d2xdt2 + 2*dxdt.*dxdt + 2*yt.*d2ydt2 + 2*dydt.*dydt));

rfac = 2*a;
theta_t = rfac*atan2(yt,xt);
dtheta_t_dt = rfac*(xt.*dydt - yt.*dxdt)./(xt.*xt + yt.*yt);
d2theta_t_dt2 = rfac*(xt.*yt.*(2*dxdt.*dxdt + yt.*d2ydt2 -2*dydt.*dydt) - ...
        xt.*xt.*(yt.*d2xdt2 + 2*dxdt.*dydt) + ...
         yt.*yt.*(2*dxdt.*dydt - yt.*d2xdt2) + xt.*xt.*xt.*d2ydt2)./ ...
     (xt.*xt + yt.*yt).^2;

xs = rho_t.*cos(theta_t);
ys = rho_t.*sin(theta_t);
dxs = drho_t_dt.*cos(theta_t) - rho_t.*sin(theta_t).*dtheta_t_dt;
dys = drho_t_dt.*sin(theta_t) + rho_t.*cos(theta_t).*dtheta_t_dt;

d2xs = d2rho_t_dt2.*cos(theta_t) - 2*drho_t_dt.*sin(theta_t).*dtheta_t_dt -...
          rho_t.*(cos(theta_t).*dtheta_t_dt.^2 + sin(theta_t).*d2theta_t_dt2);
d2ys = d2rho_t_dt2.*sin(theta_t) + 2*drho_t_dt.*cos(theta_t).*dtheta_t_dt +...
          rho_t.*(-sin(theta_t).*dtheta_t_dt.^2 + cos(theta_t).*d2theta_t_dt2);

r = [(xs(:)).'; (ys(:)).'];
d = [(dxs(:)).'; (dys(:)).'];
d2 = [(d2xs(:)).'; (d2ys(:)).'];

end
