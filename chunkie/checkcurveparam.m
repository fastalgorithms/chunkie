function dim = checkcurveparam(fcurve,ta,nout)
%CHECKCURVEPARAM check if a curve parameterization is correctly formatted 
% and return the dimension of the curve.
%
% Warning: this does *not* check that the derivatives are correct
%
% Syntax: dim = checkcurveparam(fcurve,ta)
%
% Input: 
%   fcurve - function handle of the form
%               r = fcurve(t);
%            where r is a size [dim,size(t)] arrays describing
%            the position of a curve in dim dimensions parameterized by t.
%
%            optionally, the function can be of the form 
%               [r,d] = fcurve(t);  or [r,d,d2] = fcurve(t);
%            where d is the first derivative of r with respect to t and 
%            d2 is the second derivative. 
%   ta - sample point t in parameterization to test
%   nout - number of outputs of fcurve
%
% Output:
%   dim - the dimension of the curve
%
% Examples:
%   dim = checkcurveparam(@(t) starfish(t),0.0);


if nargin < 3
    nout = 3;
end

out = cell(nout,1);
[out{:}] = fcurve(ta);

dims = zeros(nout,1);

for i = 1:nout
    sz = size(out{i});
    dims(i) = sz(1);
    assert(prod(sz(2:end)) == numel(ta), ...
        'size of each curve output should match input');
end

dim = dims(1);
assert(all(dims == dim),'dimension of curve output should be consistent');
