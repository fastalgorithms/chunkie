function ts = tangents(cgrph)
%TANGENTS return unit tangent vectors along the chunkers in the chunkgraph
%
% Syntax: ts = tangents(cgrph)
%
% Input:
%   cgrph - chunkgraph object
%
% Output:
%   ts - unit tangent vector along curve
%
% Examples:
%   ts = tangents(chnkr);
%

% author: Travis Askham (askhamwhat@gmail.com)

d = cgrph.d;

dd = sqrt(sum(abs(d).^2,1));
ts = bsxfun(@rdivide,cgrph.d,dd);

end

