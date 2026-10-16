function tup = tuples(dim,k)
%CHNK.RADPP.TUPLES index tuples of the derivative components of order k
%
% Syntax: tup = chnk.radpp.tuples(dim,k)
%
% tup is (ncomp,k), the nondecreasing tuples of 1:dim in lexicographic
% order, i.e. the order of the derivative components of a radial function
% in CHNK.RADPP.ROTATE, e.g. for dim = 2, k = 2 the rows (1,1), (1,2),
% (2,2) for G_xx, G_xy, G_yy.
%
% see also CHNK.RADPP.ROTATE

if k == 0
    tup = zeros(1,0);
else
    tup = nchoosek(1:dim+k-1,k) - (0:k-1);
end

end
