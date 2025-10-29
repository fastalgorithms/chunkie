function chnkrout = merge(chnkrs,pref)
%MERGE combine array of chunker objects into one chunker
% 
% input:
%   chnkrs - array of chunker objects, must all have same order chunks
%   pref - optional, chunkerpref object 
% output:
%   chnkrout - chunker object containing all nodes in chunker array. the
%       ordering of nodes in chnkrout has the nodes from chnkrs(1) first,
%       then chnkrs(2), etc. Adjacency information is updated as
%       appropriate to the new indices. chunker data rows are copied as
%       well. if chnkrs have different numbers of data rows, then those
%       with fewer data rows are padded with zeros on merge. 
%

if isempty(chnkrs)
  chnkrout = chunker();
  return
end
assert(isa(chnkrs,'chunker'), 'input must be of chunker type');
if nargin < 2
    pref = [];
end

% mandatory setting
pref.k = chnkrs(1).k;
t = chnkrs(1).tstor;
w = chnkrs(1).wstor;
chnkrout = chunker(pref,t,w);

for i = 1:numel(chnkrs)
  chnkrtemp = chnkrs(i);
  assert(chnkrtemp.dim == chnkrout.dim && chnkrtemp.k == chnkrout.k,...
      'chunkers to merge must be in same dimension and same order');
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


    
