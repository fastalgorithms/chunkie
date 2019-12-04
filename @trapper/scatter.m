function scatter(obj,varargin)
%SCATTER scatter plot of the coordinates of the trapper object. 
% Must be a 2D trapper
%
% see also SCATTER3

assert(obj.dim == 2,'for scatter plot must be 2D trapper');
xs = obj.r(1,:); xs = xs(:);
ys = obj.r(2,:); ys = ys(:);

scatter(xs,ys,varargin{:});