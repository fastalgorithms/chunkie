function tau = taus(chnkr)
%TAUS return unit tangent vectors along chunker curve
%
% Syntax: tau = taus(chnkr)
%
% Input:
%   chnkr - chunker object
%
% Output:
%   tau - unit tangent vector along curve
%
% Examples:
%   tau = taus(chnkr);
%

% author: Travis Askham (askhamwhat@gmail.com)

warning('taus is deprecated and will be removed. use tangents instead');

k = chnkr.k;
nch = chnkr.nch;
d = chnkr.d;

dd = sqrt(sum(abs(d).^2,1));
tau = bsxfun(@rdivide,chnkr.d,dd);

end

