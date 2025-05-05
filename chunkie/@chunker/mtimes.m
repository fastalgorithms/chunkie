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
   
end