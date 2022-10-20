function ctrs = centroids(chnkr)
%CENTROIDS
%
% Get the center of mass of each chunk

k = chnkr.k;
nch = chnkr.nch;
dim = chnkr.dim;
[~,w] = lege.exps(k);
ww = repmat( (w(:)).',dim,1);
wall = reshape(ww(:)*ones(1,nch),dim,k,nch);

ctrs = reshape(sum( (chnkr.r).*wall , 2),dim,nch)/2;
    
end