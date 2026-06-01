function qmesh = get_mesh(umesh, nchs, k)
    [x, w] = lege.exps(k);
    
    nchtot = sum(nchs(:));
    n = nchtot*k;
    qmesh.r = zeros(2,n);
    qmesh.n = zeros(2,n);
    qmesh.pseudo_normals = zeros(2,n);
    qmesh.wts = zeros(n,1);
    
    [~, ne] = size(umesh.verts);
    verts_ext = [umesh.verts, umesh.verts(:,1)];
    pseudo_normals_ext = [umesh.pseudo_normals, umesh.pseudo_normals(:,1)];
    istart = 1;
    for i=1:ne
        [xext, wext] = get_nodes(x, w, nchs(i));
        ll = length(xext);
        iend = istart+ll-1;
        qmesh.r(:,istart:iend) = verts_ext(:,i) + ...
             (verts_ext(:,i+1) - verts_ext(:,i)).*xext.';
        qmesh.n(:,istart:iend) = repmat(umesh.face_normals(:,i),[1,ll]);
        qmesh.pseudo_normals(:,istart:iend) = ...
             repmat(umesh.pseudo_normals(:,i), [1,ll]) + ...
             (pseudo_normals_ext(:,i+1) - pseudo_normals_ext(:,i)).*xext.';
        qmesh.wts(istart:iend) = wext.*umesh.lengths(i);
        istart = istart + ll;
    end

end


function [xext, wext] = get_nodes(x, w, nch)
    tabs = 0:1/nch:1;
    tas = tabs(1:nch);
    tbs = tabs(2:end);
    xext = (x(:)+1)/2.*(tbs - tas) + tas;
    xext = xext(:);
    wext = w(:)/2.*(tbs-tas);
    wext = wext(:);    
end