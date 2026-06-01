function h = minus(f,g)
% - pointwise subtraction for chunkgraph class
%
% Currently only supported for subtracting a constant vector from the
% chunkgraph, translating the entire object

if ~isa(f,'chunkgraph')
    % vector - chunkgraph = -(chunkgraph - vector)
    h = -1*minus(g,f);
    return
elseif isa(g,'numeric')
    h = f;
    for j = 1:numel(h.echnks)
        h.echnks(j) = h.echnks(j) - g(:);
    end
    if (size(h.verts,2) < 1)
        return
    end
    assert(numel(g) == size(h.verts,1),...
        'CHUNKGRAPH:minus vector must have same dimension as chunkgraph points');
    h.verts = h.verts - g(:);
else
    error('CHUNKGRAPH:minus:invalid', ...
       'minus only defined for a chunkgraph object minus a vector of appropriate size');
end
