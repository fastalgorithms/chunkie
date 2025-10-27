function cgrph = tochunkgraph(chnkr)
% TOCHUNKGRAPH converts a chunker into a chunkgraph with one edge per
% connected component

chnkr = sort(chnkr);
[~,~,info] = sortinfo(chnkr);

ncomp = info.ncomp;
nchs = info.nchs;

istart = 1;
nverts = 0;
verts = [];
edge2verts = [];
cparams = cell(1,ncomp);
fchnks = cell(1,ncomp);

for i = 1:ncomp
    nch = nchs(i);
    vs = chunkends(chnkr,[istart,istart+nch-1]);
    
    ifclosed = info.ifclosed(i);
    if ifclosed
        vs = vs(:,1);
        verts = [verts, vs];
        edge2verts = [edge2verts, [nverts+1;nverts+1]];
        nverts = nverts + 1;
    else
        vs = vs(:, [1,4]);
        verts = [verts, vs];
        edge2verts = [edge2verts, [nverts+1;nverts+2]];
        nverts = nverts + 2;
    end

    rchnkr = []; 
    rchnkr.r  = chnkr.r(:,:,istart:(istart+nch-1));
    rchnkr.d  = chnkr.d(:,:,istart:(istart+nch-1));
    rchnkr.d2 = chnkr.d2(:,:,istart:(istart+nch-1));
    fchnks{i} = chunkerpoints(rchnkr,struct('ifclosed',ifclosed));
    cparams{i}.ifclosed = ifclosed;
    istart = istart + nch;
end

pref = [];

pref.dim = chnkr.dim;
pref.k   = chnkr.k;
cgrph = chunkgraph(verts,edge2verts,fchnks,cparams,pref);

end