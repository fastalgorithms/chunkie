
function [proxy,pnorm,pw] = proxy_circ_pts(p)

if nargin < 1
    p = 64;
end

theta = (0:(p-1))*2*pi/p;

proxy = [1.5*cos(theta); 1.5*sin(theta)];

pnorm = [cos(theta); sin(theta)];

pw = ones(p,1)*(2.0*pi*1.5)/p;

end