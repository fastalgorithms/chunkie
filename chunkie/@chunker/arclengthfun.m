function [s,nchs,chnkr] = arclengthfun(chnkr)
%ARCLENGTHFUN
% returns the values of s (arclength) along the
% curve. if curve is multiple components, get
% a different function per component
%
  
[~,~,info] = sortinfo(chnkr);
chnkr = sort(chnkr);

ncomp = info.ncomp;
nchs = info.nchs;
istart = 1;

A = lege.intmat(chnkr.k);
ds = arclengthdens(chnkr);
chunklens = chunklen(chnkr);
s = A*ds;

for i = 1:ncomp
    nch = nchs(i);
    indsi = istart:(istart+nch-1);
    sstart = 0;
    for j = indsi
        lj = chunklens(j);
        s(:,j) = s(:,j)+sstart;
        sstart = sstart+lj;
    end
    istart = istart + nch;
end
