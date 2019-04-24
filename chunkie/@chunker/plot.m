function plot(obj,varargin)
%PLOT plot the coordinates of the chunker object
% Must be a 2D chunker

assert(obj.dim == 2,'for plot must be 2D chunker');

[obj,ifclosed] = sort(obj);
xs = obj.r(1,:,:); xs = xs(:);
ys = obj.r(2,:,:); ys = ys(:);

if ifclosed
    xs = [xs(:); xs(1)];
    ys = [ys(:); ys(1)];
end

plot(xs,ys,varargin{:})