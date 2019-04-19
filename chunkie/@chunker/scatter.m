function scatter(obj,varargin)
%SCATTER scatter plot of the coordinates of the chunker object. 
% Must be a 2D chunker
%
% see also SCATTER3

assert(obj.dim == 2,'for scatter plot must be 2D chunker');
xs = obj.r(1,:,:); xs = xs(:);
ys = obj.r(2,:,:); ys = ys(:);

scatter(xs,ys,varargin{:});