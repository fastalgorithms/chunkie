function [pr,ptau,pw,pin] = proxy_rect_pts(lxy, npxy, opts)
%PROXY_SQUARE_PTS return the proxy points, unit tangents, weights, and 
% function handle for determining if points are within proxy surface.
% This function is for the proxy surface around a unit box centered at 
% the origin. The proxy points lie on the 
% [-lxy(1),lxy(1)]*[-lxy(2), lxy(2)] rectangle
%
% Input: 
%   lxy(2) - half edge lengths in x and y direction
%   npxy(2) - number of proxy points in x and y direction
%   opts - optional input
%     opts.iflege (false) whether to use legendre or equispaced nodes
%   
% Output:
%   pr - (2,2*(npxy(1) + npxy(2))) array of proxy point locations 
%          


if nargin < 1
    lxy = [1,1];
end

if nargin < 2
    npxy = [10, 10];
end

if nargin < 3
    opts = [];
end

iflege = 0;
if isfield(opts, 'iflege')
    iflege = opts.iflege;
end

nx = npxy(1);
ny = npxy(2);

lx = lxy(1);
ly = lxy(2);


if ~iflege
    pts_x = -lx +2*lx*(0:(nx-1))*1.0/nx;
    pw_x = ones(nx,1)*(2*lx/nx);

    pts_y = -ly +2*ly*(0:(ny-1))*1.0/ny;
    pw_y = ones(ny,1)*(2*ly/ny);
else
    [pts_x, pw_x] = lege.exps(nx);
    [pts_y, pw_y] = lege.exps(ny);
    pts_x = pts_x.'*lx;
    pw_x = pw_x*lx;

    pts_y = pts_y.'*ly;
    pw_y = pw_y*ly;

end


one_x = ones(size(pts_x));
one_y = ones(size(pts_y));

pr = [pts_x, one_y*lx, -pts_x, -lx*one_y;
    -ly*one_x, pts_y, ly*one_x, -pts_y];

pw = [pw_x; pw_y; pw_x; pw_y];

ptau = [one_x, one_y*0, -one_x, one_y*0;
    one_x*0, one_y, one_x*0, -one_y];

pin = @(x) max(abs(x./lpxy),[],1) < 1;

end