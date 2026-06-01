function g = plus(f,g)
% + pointwise addition for chunkgraph class
% 
% Currently only supported for adding a constant vector to the chunkgraph,
% translating the entire object

if ~isa(g,'chunkgraph')
    g = plus(g,f);
elseif  isa(f,'numeric')
    for j = 1:numel(g.echnks)
        g.echnks(j) = f(:) + g.echnks(j);
    end
    if (size(g.verts,2) < 1)
        return
    end
    assert(numel(f) == size(g.verts,1),...
        'CHUNKGRAPH:plus vector must have same dimension as chunkgraph points');
    g.verts = f(:) + g.verts;
else
    error('CHUNKGRAPH:plus:invalid', ...
       'plus only defined for sum of an appropriate size vector and chunkgraph object');
end