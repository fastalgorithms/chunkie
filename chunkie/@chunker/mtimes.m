function chnkr = mtimes(A,chnkr)
% * Matrix multiplication for chunker. Transforms curve positions by 
% r -> A*r and updates derivatives and normals accordingly. Note that this
% routine can change the orientation if det(A) is not positive.
%

if isa(A,"numeric") && isa(chnkr,"chunker")
    dim = size(chnkr.r,1);
    if isscalar(A)
        A = A*eye(dim);
    end
    [m,n] = size(A);
    if n == dim && m == dim
        chnkr.r(:,:) = A*chnkr.r(:,:);
        chnkr.d(:,:) = A*chnkr.d(:,:);
        chnkr.d2(:,:) = A*chnkr.d2(:,:);
        chnkr.n = normals(chnkr);
        chnkr.wts = weights(chnkr);
    else
        error("CHUNKER:mtimes:invalid",...
            "matrix must have compatible size for transforming coordinates");
    end
elseif isnumeric(chnkr) && isscalar(chnkr)
    chnkr = mtimes(chnkr,A);
elseif isnumeric(chnkr)
    error("CHUNKER:mtimes:invalid",...
            "product of chunker and matrix only defined for matrix on left");
else
    error("CHUNKER:mtimes:invalid",...
        "type not supported for product with chunker");
end