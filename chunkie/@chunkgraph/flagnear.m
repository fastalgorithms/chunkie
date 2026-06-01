function flag = flagnear(obj,pts,opts)
%FLAGNEAR flag points which are within a chunklength (or given factor of
% a chunklength) of any point on each chunk. on return is a sparse
% logical array. If the (i,j) entry is non zero then pts(:,i) is close
% to chunk j in chunkgraph.
%
% Syntax: flag = flagnear(obj,pts,opts)
%
% Input:
%   obj - chunkgraph object describing curve
%   pts - (obj.dim,nt) array of point coordinates
%
% Optional input:
%   opts - options structure
%       opts.fac = factor of chunklength to check distance against (1.0)
%
% Output:
%   flag - (nt,nch) sparse array. a non zero entry (i,j) means 
%       that the distance from pts(:,i) to at least one node on 
%       chunk j is less than opts.fac*length of chunk j.
%   where nch is thte total number of chunks int he chunkgraph
%

% author: Travis Askham (askhawhat@gmail.com)

chnkrtotal = merge(obj.echnks);
flag = chnkrtotal.flagnear(pts, opts);

end