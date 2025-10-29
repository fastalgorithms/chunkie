function [cgrph, isort] = split_chunkgraph(cgrph, iedge, chsplit)

chnkr_loc = cgrph.echnks(iedge); nch = chnkr_loc.nch;
npts = [cgrph.echnks.npt];

if (chsplit == 0) || (chsplit == nch)
    isort = 1:sum(npts);
    return
end

opts = []; opts.ifclosed = 0;
src1 = [];
src1.r = chnkr_loc.r(:,:,1:chsplit);
src1.d = chnkr_loc.d(:,:,1:chsplit);
src1.d2 = chnkr_loc.d2(:,:,1:chsplit);

src2 = [];
src2.r = chnkr_loc.r(:,:,(chsplit+1):nch);
src2.d = chnkr_loc.d(:,:,(chsplit+1):nch);
src2.d2 = chnkr_loc.d2(:,:,(chsplit+1):nch);

chnkr1 = chunkerpoints(src1, opts);
chnkr2 = chunkerpoints(src2, opts);

vert_new = chunkends(chnkr1, chsplit); vert_new = vert_new(:,2);

nverts = size(cgrph.verts,2);
verts = [cgrph.verts, vert_new];

% get the edges
edgesendverts_new = cgrph.edgesendverts;

iverts = edgesendverts_new(:,iedge);

edgesendverts_new(:,iedge) = [iverts(1); nverts+1];
edgesendverts_new = [edgesendverts_new, [nverts+1; iverts(2)]];

echnks = cgrph.echnks;
echnks(iedge) = chnkr1;
echnks(size(echnks,2)+1) = chnkr2;


% assemble the new graph
cgrph.verts = verts;
cgrph.edgesendverts = edgesendverts_new;
cgrph.v2emat        = build_v2emat(cgrph);
cgrph.echnks        = echnks;
cgrph.vstruc = procverts(cgrph);
cgrph.wts = weights(cgrph);

cgrph.regions = findregions(cgrph);

% get sorting info
npt_tot = cgrph.npt; npt2 = chnkr2.npt;
isortl = 1:sum([npts(1:iedge-1), chnkr1.npt]);
isortm = (npt_tot-npt2+1):npt_tot;
isortr = isortl(end)+1:(npt_tot-npt2);

itmp = [isortl, isortm, isortr];
isort = 1:npt_tot;
isort(itmp) =  1:npt_tot;
end