function dim = checkcurveparam(fcurve,ta)
%CHECKCURVEPARAM
%

nout = 3;
out = cell(nout,1);
[out{:}] = fcurve(ta);

dims = zeros(nout,1);

for i = 1:nout
    sz = size(out{i});
    dims(i) = sz(1);
    assert(prod(sz(2:end)) == numel(ta),'size of each curve output should match input');
end

dim = dims(1);
assert(all(dims == dim),'dimension of curve output should be consistent');
