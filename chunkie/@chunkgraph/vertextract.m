function [inds,isgn] = vertextract(ivert,cgrph)

irights = find(cgrph.edge2verts(:,ivert) ==  1);
ilefts  = find(cgrph.edge2verts(:,ivert) == -1);

drights = [];
for i=1:numel(irights)
drights = [drights,-cgrph.echnks(irights(i)).d(:,end,end)];
end

dlefts = [];
for i=1:numel(ilefts)
dlefts = [dlefts,cgrph.echnks(ilefts(i)).d(:,1,1)];
end

ds = [drights,dlefts];
angs = atan2(ds(2,:),ds(1,:));

[angs,isrtededges] = sort(angs);
inds = [irights',ilefts'];
isgn = [ones(size(irights')),-ones(size(ilefts'))];
inds = inds(isrtededges);
isgn = isgn(isrtededges);
end
