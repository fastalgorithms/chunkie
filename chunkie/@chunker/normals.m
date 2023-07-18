function rnorms = normals(chnkr)
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

assert(chnkr.dim == 2,'normals only implemented for dim=2');
k = chnkr.k;
nch = chnkr.nch;

dd = sqrt(chnkr.dstor(1,:,1:nch).^2 + chnkr.dstor(2,:,1:nch).^2);
rnorms = zeros(2,k,nch);
rnorms(1,:,:) = chnkr.dstor(2,:,1:nch)./dd;
rnorms(2,:,:) = -chnkr.dstor(1,:,1:nch)./dd;

end

