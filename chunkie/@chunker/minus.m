function h = minus(f,g)
% - pointwise subtraction for chunker class
%
% Currently only supported for subtracting a constant vector from the chunker,
% translating the curve

if ~isa(f,'chunker')
    % vector - chunker = -(chunker - vector)
    h = -1*minus(g,f);
elseif isa(g,'numeric') && numel(g) == size(f.r,1)
    h = f;
    h.r(:,:) = f.r(:,:) - g(:);
else
    error('CHUNKER:minus:invalid', ...
       'minus only supported for a chunker object minus a vector of appropriate size');
end
