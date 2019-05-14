function [pr,ptau,pw,pin] = proxy_square_pts(porder)
%PROXY_SQUARE_PTS return the proxy points, unit tangents, weights, and 
% function handle for determining if points are within proxy surface.
% This function is for the proxy surface around a unit box centered at 
% the origin. The proxy points lie on the [-1.5,1.5]^2 square
%
% Input: 
%   porder - number of points on proxy surface
%   
% Output:
%   pr - (2,porder) array of proxy point locations (centered at 0 for a 
%           box of side length 1)


if nargin < 1
    porder = 64;
end

assert(mod(porder,4) == 0, ...
    'number of proxy points on square should be multiple of 4')

po4 = porder/4;

pts = -1.5+3*(0:(po4-1))*1.0/po4;
one = ones(size(pts));

pr = [pts, one*1.5, -pts, -1.5*one;
    -1.5*one, pts, 1.5*one, -pts];

ptau = [one, one*0, -one, one*0;
    one*0, one, one*0, -one];

pw = ones(porder,1)*(4.0*3.0/porder);

pin = @(x) max(abs(x),[],1) < 1.5;

end