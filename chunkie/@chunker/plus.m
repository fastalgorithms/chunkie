function chnkr = plus(v,chnkr)
% + pointwise addition for chunker class
% 
% Currently only supported for adding a constant vector to the chunker,
% translating the curve 

if isa(v,'numeric') && isa(chnkr,'chunker')
    assert(numel(v) == size(chnkr.r,1),'CHUNKER:plus sizes incompatible');
    chnkr.r(:,:) = chnkr.r(:,:) + v(:);
else
    error('CHUNKER:plus:invalid', ...
       'v must be a vector of appropriate size and chnkr a chunker object');
end