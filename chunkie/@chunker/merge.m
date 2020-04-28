function chnkrout = merge(chnkrs)

if isempty(chnkrs)
    chnkrout = chunker();
    return
end
assert(isa(chnkrs,'chunker'), 'input must be of chunker type');

chnkrout = chnkrs(1);

for i = 2:length(chnkrs)
    chnkrtemp = chnkrs(i);
    assert(chnkrtemp.dim == chnkrout.dim,...
        'chunkers to merge must be in same dimension');
    nch = chnkrtemp.nch;
    nchold = chnkrout.nch;
    istart = nchold+1;
    iend = istart+nch-1;
    
    chnkrout = chnkrout.addchunk(nch);
    chnkrout.r(:,:,istart:iend) = chnkrtemp.r;
    chnkrout.d(:,:,istart:iend) = chnkrtemp.d;
    chnkrout.d2(:,:,istart:iend) = chnkrtemp.d2;
    chnkrout.adj(:,istart:iend) = chnkrtemp.adj + nchold;
    chnkrout.h(istart:iend) = chnkrtemp.h;
    
end

    
