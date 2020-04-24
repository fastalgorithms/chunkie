function ts = tangents(chnkr)
%TANGENTS return unit tangent vectors along chunker curve
%
% Syntax: ts = tangents(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   ts - unit tangent vector along curve
%
% Examples:
%   ts = tangents(chnkr);
%

% author: Travis Askham (askhamwhat@gmail.com)

d = chnkr.d;

dd = sqrt(sum(abs(d).^2,1));
ts = bsxfun(@rdivide,chnkr.d,dd);

end

