function flag = flagnear_rectangle_grid(obj,x,y,opts)
%FLAGNEAR_RECTANGLE_GRID flag points which require special quadrature
% (or given factor of a chunklength) of any point on each chunk. 
% On return is a sparse logical array. If the (i,j) entry is non zero 
% then pts(:,i) is close to chunk j in chnkr.
%
% Syntax: flag = flagnear(obj,pts,opts)
%
% Input:
%   obj - chunkgraph object describing curve
%   x - nx array of x coordinates for tensor product grid
%   y - ny array of y coordinates for tensor product grid
%
% Grid points are interpreted in the order returned by 
% MATLAB's meshgrid function [xx,yy] = meshgrid(x,y).
%
% Optional input:
%   opts - options structure
%       opts.rho = Bernstein ellipse parameter (default=1.8)
%
% Output:
%   flag - (nx*ny, nch) sparse array. a non zero entry (i,j) means 
%       that the distance from pts(:,i) to at least one node on 
%       chunk j is less than opts.fac*length of chunk j.
% 
%   Here nch is the total number of chunks on the chunkgraph
%

% author: Travis Askham (askhawhat@gmail.com)

chnkrtotal = merge(obj.echnks);
flag = chnkrtotal.flagnear_rectangle_grid(x, y, opts);

end