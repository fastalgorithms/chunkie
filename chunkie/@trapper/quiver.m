function quiver(obj,varargin)
%QUIVER quiver plot of chunker normal vectors in 2 dimensions

assert(obj.dim == 2,'for quiver plot must be 2D chunker');

xs = obj.r(1,:); xs = xs(:);
ys = obj.r(2,:); ys = ys(:);

rnorms = normals(obj);

u = rnorms(1,:); u = u(:);
v = rnorms(2,:); v = v(:);

quiver(xs,ys,u,v,varargin{:});