function [umesh,dmesh,qmesh,scales,levels] = merge(umeshes,dmeshes,...
    qmeshes,src_codes,targ_codes)

    nv = 0;
    for ii=1:numel(umeshes)
        nv = nv + size(umeshes{ii}.verts,2);
    end

    verts = zeros(2,nv);
    centroids = zeros(2,nv);
    lengths = zeros(1,nv);
    face_normals = zeros(2,nv);
    pseudo_normals = zeros(2,nv);

    ind = 1;
    for ii=1:numel(umeshes)
        nv_loc = size(umeshes{ii}.verts,2);
        ind2 = ind + nv_loc - 1;
        verts(:,ind:ind2)  = umeshes{ii}.verts;
        centroids(:,ind:ind2)  = umeshes{ii}.centroids;
        lengths(:,ind:ind2) = umeshes{ii}.lengths;
        face_normals(:,ind:ind2)  = umeshes{ii}.face_normals;
        pseudo_normals(:,ind:ind2)  = umeshes{ii}.pseudo_normals;
        ind = ind2 + 1;
    end

    umesh = [];
    umesh.verts = verts;
    umesh.centroids = centroids;
    umesh.lengths = lengths;
    umesh.face_normals = face_normals;
    umesh.pseudo_normals = pseudo_normals;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    npd = 0;
    for ii=1:numel(dmeshes)
        npd = npd + size(dmeshes{ii}.r,2);
    end

    r = zeros(2,npd);
    n = zeros(2,npd);
    pseudo_normals = zeros(2,npd);
    wts = zeros(npd,1);
    levels = ones(npd,1);
    
    ind = 1;
    for ii=1:numel(dmeshes)
        np_loc = size(dmeshes{ii}.r,2);
        ind2 = ind + np_loc - 1;
        r(:,ind:ind2) = dmeshes{ii}.r;
        n(:,ind:ind2) = dmeshes{ii}.n;
        pseudo_normals(:,ind:ind2) = dmeshes{ii}.pseudo_normals;
        wts(ind:ind2,1) = dmeshes{ii}.wts;
        if (nargin>4)
            levels(ind:ind2) = targ_codes(ii);
        end
        ind = ind2 + 1;
    end

    dmesh = [];
    dmesh.r = r;
    dmesh.n = n;
    dmesh.pseudo_normals = pseudo_normals;
    dmesh.wts = wts;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nqd = 0;
    for ii=1:numel(qmeshes)
        nqd = nqd + size(qmeshes{ii}.r,2);
    end

    r = zeros(2,nqd);
    n = zeros(2,nqd);
    pseudo_normals = zeros(2,nqd);
    wts = zeros(nqd,1);
    scales = ones(nqd,1);
    
    ind = 1;
    for ii=1:numel(qmeshes)
        np_loc = size(qmeshes{ii}.r,2);
        ind2 = ind + np_loc - 1;
        r(:,ind:ind2) = qmeshes{ii}.r;
        n(:,ind:ind2) = qmeshes{ii}.n;
        pseudo_normals(:,ind:ind2) = qmeshes{ii}.pseudo_normals;
        wts(ind:ind2,1) = qmeshes{ii}.wts;
        if (nargin>3)
            scales(ind:ind2) = src_codes(ii);
        end
        ind = ind2 + 1;
    end

    qmesh = [];
    qmesh.r = r;
    qmesh.n = n;
    qmesh.pseudo_normals = pseudo_normals;
    qmesh.wts = wts;

end