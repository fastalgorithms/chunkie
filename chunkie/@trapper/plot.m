function plot(obj,varargin)
%PLOT plot the coordinates of the trapper object
% Must be a 2D trapper

assert(obj.dim == 2,'for plot must be 2D trapper');

xs = obj.r(1,:); ys = obj.r(2,:);
xs = [xs, xs(1)]; ys = [ys, ys(1)];
plot(xs,ys,varargin{:})