function g = plus(f,g)
% + pointwise addition for chunker class
% 
% Currently only supported for adding a constant vector to the chunker,
% translating the curve 

if ~isa(g,"chunker")
    g = plus(g,f);
elseif isa(f,'numeric') && numel(f) == size(g.r,1)
    assert(numel(f) == size(g.r,1),'CHUNKER:plus sizes incompatible');
    g.r(:,:) = g.r(:,:) + f(:);
else
    error('CHUNKER:plus:invalid', ...
       'plus only supported for a vector of appropriate size and a chunker object');
end