function [umesh] = get_umesh(verts)
    umesh = [];
    umesh.verts = verts;
    [~, nv] = size(verts);
    verts_ext = [verts, verts(:,1)];
    d = verts_ext(:,2:end) - verts_ext(:,1:nv);
    umesh.centroids = 0.5*(verts_ext(:,2:end) + verts_ext(:,1:nv));
    dd = sqrt(sum(d.^2, 1));
    umesh.lengths = dd;
    face_normals = zeros(2,nv);
    face_normals(1,:) = d(2,:)./dd;
    face_normals(2,:) = -d(1,:)./dd;

    umesh.face_normals = face_normals;
    face_normals_ext = [face_normals(:,end), face_normals];
    pseudo_normals = 0.5*(face_normals_ext(:,2:end) + ...
               face_normals_ext(:,1:nv));
    dd = sqrt(sum(pseudo_normals.^2, 1));
    pseudo_normals = pseudo_normals./dd;

    umesh.pseudo_normals = pseudo_normals;
    
end