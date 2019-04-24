function a = area(chnkr)
%AREA compute area of a closed, 2D chunker curve

assert(chnkr.dim == 2,'area only well-defined for 2d curves');
assert(nnz(chnkr.adj < 1) == 0, 'area only well-defined for closed curves');

wts = whts(chnkr);
rnorm = normals(chnkr);
a = sum(sum(wts.*squeeze(sum(rnorm.*(chnkr.r),1))))/chnkr.dim;

end