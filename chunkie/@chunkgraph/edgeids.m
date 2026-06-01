function ids = edgeids(cgrph, ichs)
%edgeids return the indices corresponding to points in cgrph.echnks(ichs)
npts = [cgrph.echnks.npt];
% get starting index for each edge
idstart = [1,cumsum(npts)+1];

ids = [];
for i = 1:length(ichs)
    ids = [ids, idstart(ichs(i)):(idstart(ichs(i)+1)-1)];
end
end