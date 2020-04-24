function [y] = planewave(kvec,r)

y=exp(1i*sum(bsxfun(@times,kvec(:),r(:,:))));