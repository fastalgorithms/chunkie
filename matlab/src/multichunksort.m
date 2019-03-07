function chunkerout = multichunksort(chunker,nchs,ifexplode)
%MULTICHUNKSORT takes in a description of multiple chunked 
% pieces of domain (with the number of chunks on each in nchs)
% and sorts each piece. Useful for plotting
%
% if ifexplode == true, then the output is an array with
% a chunker for each piece of domain (default false)

checkchunker(chunker);

if nargin < 3
    ifexplode = false;
end

nch = chunker.nch;
k = chunker.k;
chunker.chunks = reshape(chunker.chunks,2,k,nch);
chunker.ders= reshape(chunker.ders,2,k,nch);
chunker.ders2= reshape(chunker.ders2,2,k,nch);
chunker.adjs= reshape(chunker.adjs,2,nch);
chunker.hs= reshape(chunker.hs,nch,1);

if ifexplode
    chunkerout = cell(length(nchs),1);
    ind = 1;
    for i = 1:length(nchs)
        chunkerout{i} = [];
        nchi = nchs(i);
        chunkerout{i}.k = chunker.k;
        chunkerout{i}.nch = nchi;
        chunkerout{i}.chunks = zeros(2,k,nchi);
        chunkerout{i}.ders = zeros(2,k,nchi);
        chunkerout{i}.ders2 = zeros(2,k,nchi);
        chunkerout{i}.adjs = zeros(2,nchi);
        chunkerout{i}.hs = zeros(nchi,1);
        chunkerout{i}.chunks(:,:,:) = chunker.chunks(:,:,ind:ind+nchi-1);
        chunkerout{i}.ders(:,:,:) = chunker.ders(:,:,ind:ind+nchi-1);
        chunkerout{i}.ders2(:,:,:) = chunker.ders2(:,:,ind:ind+nchi-1);
        chunkerout{i}.adjs(:,:) = chunker.adjs(:,ind:ind+nchi-1)-ind+1;
        chunkerout{i}.hs(:) = chunker.hs(ind:ind+nchi-1);
        checkchunker(chunkerout{i});
        chunkerout{i} = chunksort(chunkerout{i});
        ind = ind+nchi;
    end
    
else
    chunkerout = [];
    chunkerout.chunks = chunker.chunks;
    chunkerout.ders = chunker.ders;
    chunkerout.ders2 = chunker.ders2;
    chunkerout.adjs = chunker.adjs;
    chunkerout.hs = chunker.hs;
    ind = 1;
    for i = 1:length(nchs)
        chunkert = [];
        nchi = nchs(i);
        chunkert.k = chunker.k;
        chunkert.nch = nchi;
        chunkert.chunks = chunker.chunks(:,:,ind:ind+nchi-1);
        chunkert.ders = chunker.ders(:,:,ind:ind+nchi-1);
        chunkert.ders2 = chunker.ders2(:,:,ind:ind+nchi-1);
        chunkert.adjs = chunker.adjs(:,ind:ind+nchi-1);
        chunkert.hs = chunker.hs(ind:ind+nchi-1);
        chunkert = chunksort(chunkert);
        chunkerout.chunks(:,:,ind:ind+nchi-1) = chunkert.chunks;
        chunkerout.ders(:,:,ind:ind+nchi-1) = chunkert.ders;
        chunkerout.ders2(:,:,ind:ind+nchi-1) = chunkert.ders2 ;
        chunkerout.adjs(:,ind:ind+nchi-1) = chunkert.adjs;
        chunkerout.hs(ind:ind+nchi-1) = chunkert.hs;
        ind = ind+nchi;
    end
    
end
        
end