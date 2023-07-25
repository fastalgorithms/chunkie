function rnorms = normals(chnkgrph)
%NORMALS compute normal vectors along boundary of 2D chunkgraph
% 
% Syntax: rnorms = normals(cgrph)
%
% Input:
%   cgrph - chunkgraph object
%
% Output:
%   rnorms - normals along each edge of the cunkgraph
%
% Examples:
%   rnorms = normals(cgrph)
%

% author: Jeremy Hoskins

% merge the edge chunker objects into a single chunk object
chnkr = merge(chnkgrph.echnks);
assert(chnkr.dim == 2,'normals only implemented for dim=2');
k = chnkr.k;
nch = chnkr.nch;
d = chnkr.d;

% compute the normals from the tangent information
dd = sqrt(d(1,:,:).^2 + d(2,:,:).^2);
rnorms = zeros(2,k,nch);
rnorms(1,:,:) = d(2,:,:)./dd;
rnorms(2,:,:) = -d(1,:,:)./dd;

end

