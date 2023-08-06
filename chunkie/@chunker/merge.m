function chnkrout = merge(chnkrs)

if isempty(chnkrs)
  chnkrout = chunker();
  return
end
assert(isa(chnkrs,'chunker'), 'input must be of chunker type');

chnkrout = chunker();
%chnkrout = chnkrs(1);


for i = 1:length(chnkrs)
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
  ipos = chnkrtemp.adj > 0;
  chnkrtemp.adj(ipos) = chnkrtemp.adj(ipos)+nchold;
  chnkrout.adj(:,istart:iend) = chnkrtemp.adj;
  chnkrout.h(istart:iend) = chnkrtemp.h;
  chnkrout.wts(:,istart:iend) = chnkrtemp.wts;
  chnkrout.n(:,:,istart:iend) = chnkrtemp.n;
end


mmax = 0;
nchold = 0;
for i=1:length(chnkrs)
    if(chnkrs(i).hasdata)
        [m,~,~] = size(chnkrs(i).data);
        if(m>mmax) 
            mmax = m;
        end
    end
end
if(mmax>0)
    chnkrout = chnkrout.makedatarows(mmax);
    for i=1:length(chnkrs)
      chnkrtemp = chnkrs(i);
  
      nch = chnkrtemp.nch;
      istart = nchold+1;
      iend = istart+nch-1;
      
      if(chnkrtemp.hasdata)
        [m,~,~] = size(chnkrtemp.data);
        chnkrout.data(1:m,:,istart:iend) = chnkrtemp.data;
      end
      nchold = nchold + nch;
    end
end


    
