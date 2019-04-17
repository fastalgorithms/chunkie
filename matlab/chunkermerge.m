function [chunker,nchs] = chunkermerge(chunkers,ifrev,ifsort)
%CHUNKERMERGE takes in a cell array of chunker objects and 
% merges them into one big chunker with nchs array indicating
% how many chunks correspond to each

assert(iscell(chunkers), ...
    'input must be cell array of chunker objs');

chunker = [];
nchs = [];

if isempty(chunkers)
    return;
end

if nargin < 2
    ifrev = false(length(chunkers),1);
end
if nargin < 3
    ifsort = false(length(chunkers),1);
end

if length(ifrev) == 1
    ifrev = ifrev(1)*true(length(chunkers),1);
end
if length(ifsort) == 1
    ifsort = ifsort(1)*true(length(chunkers),1);
end

nchs = zeros(length(chunkers),1);
nchs(1) = chunkers{1}.nch;
k = chunkers{1}.k;
for i = 2:length(chunkers)
    assert(chunkers{i}.k == k, ...
        'when merging, chunkers must have same order');
    nchs(i) = chunkers{i}.nch;
end

nch = sum(nchs);
chunks = zeros(2,k,nch); ders = zeros(2,k,nch);
ders2 = zeros(2,k,nch); adjs = zeros(2,nch); hs = zeros(nch,1);

nchtot = 0;
for i = 1:length(chunkers)
    chunkt = chunkers{i};
    if ifrev(i)
        chunkt = chunkreverse(chunkt);
    end
    if ifsort(i)
        chunkt = chunksort(chunkt);
    end
    nchi = chunkt.nch;
    chunks(:,:,nchtot+1:nchtot+nchi) = reshape(chunkt.chunks,2,k,nchi);
    ders(:,:,nchtot+1:nchtot+nchi) = reshape(chunkt.ders,2,k,nchi);
    ders2(:,:,nchtot+1:nchtot+nchi) = reshape(chunkt.ders2,2,k,nchi);
    adjs(:,nchtot+1:nchtot+nchi) = nchtot+chunkt.adjs;
    hs(nchtot+1:nchtot+nchi) = chunkt.hs;
    nchtot = nchtot + nchi;
end

chunker.k = k; chunker.nch = nch;
chunker.chunks = chunks;
chunker.adjs = adjs;
chunker.ders = ders;
chunker.ders2 = ders2;
chunker.hs = hs;  

end

    
    