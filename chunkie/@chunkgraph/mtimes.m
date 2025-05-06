function cg = mtimes(A,cg)
% * Matrix multiplication for chunkgraph. Transforms curve positions by 
% r -> A*r and updates derivatives, normals and graph vertices 
% accordingly. Note that this routine can change the orientation if 
% det(A) is not positive.

if isa(A,"numeric") && isa(cg,"chunkgraph")
    for j = 1:length(cg.echnks)
        cg.echnks(j) = A*cg.echnks(j);
    end
    
    if size(cg.verts,2) < 1
        return
    end
    dim = size(cg.verts,1);
    if isscalar(A)
        A = A*eye(dim);
    end
    [m,n] = size(A);
    if n == dim && m == dim
        cg.verts = A*cg.verts;
    else
        error("CHUNKER:mtimes:invalid",...
            "matrix must have compatible size for transforming coordinates");
    end
elseif isnumeric(cg) && isscalar(cg)
    cg = mtimes(cg,A);
elseif isnumeric(cg)
    error("CHUNKER:mtimes:invalid",...
        "product of chunkgraph and matrix only defined for matrix on left");
else
    error("CHUNKER:mtimes:invalid",...
        "type not supported for product with chunkgraph");
end