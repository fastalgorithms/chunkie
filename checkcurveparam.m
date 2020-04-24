function dim = checkcurveparam(fcurve,ta)
%CHECKCURVEPARAM check if a curve parameterization is correctly formatted 
% and return the dimension of the curve.
%
% Warning: this does *not* check that the derivatives are correct
%
% Syntax: dim = checkcurveparam(fcurve,ta)
%
% Input: 
%   fcurve - function handle of the form
%               [r,d,d2] = fcurve(t)
%            where r, d, d2 are size [dim,size(t)] arrays describing
%            position, first derivative, and second derivative of a curve
%            in dim dimensions parameterized by t.
%   ta - sample point t in parameterization to test
%
% Output:
%   dim - the dimension of the curve
%
% Examples:
%   dim = checkcurveparam(@(t) starfish(t),0.0);


nout = 3;
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
