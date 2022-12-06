function rnorms = normals(chnkgrph)
%NORMALS compute normal vectors along boundary of 2D curve
% 
% Syntax: rnorms = normals(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   rnorms - normals along boundary of chunker curve
%
% Examples:
%   rnorms = normals(chnkr)
%

% author: Travis Askham (askhamwhat@gmail.com)

chnkr = merge(chnkgrph.echnks);
assert(chnkr.dim == 2,'normals only implemented for dim=2');
k = chnkr.k;
nch = chnkr.nch;
d = chnkr.d;

dd = sqrt(d(1,:,:).^2 + d(2,:,:).^2);
rnorms = zeros(2,k,nch);
rnorms(1,:,:) = d(2,:,:)./dd;
rnorms(2,:,:) = -d(1,:,:)./dd;

end

