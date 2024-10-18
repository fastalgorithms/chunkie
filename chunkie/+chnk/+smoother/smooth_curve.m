nv = 3;
z = exp(1j*2*pi*(1:nv)/nv);
verts = [real(z); imag(z)];
nchs = randi([1,10], 1, nv);
k = 16;

kquad = 32;



umesh = chnk.smoother.get_umesh(verts);
dmesh = chnk.smoother.get_mesh(umesh, nchs, k);
qmesh = get_mesh(umesh, nchs, kquad);

figure(1)
clf
plot(verts(1,:), verts(2,:), 'k.','MarkerSize',20); hold on;
plot(umesh.centroids(1,:), umesh.centroids(2,:), 'r.', 'MarkerSize',20);
quiver(verts(1,:), verts(2,:), umesh.pseudo_normals(1,:), umesh.pseudo_normals(2,:),'k')
quiver(umesh.centroids(1,:), umesh.centroids(2,:), umesh.face_normals(1,:), umesh.face_normals(2,:),'r')

figure(2)
clf
plot(dmesh.r(1,:), dmesh.r(2,:),'k.'); hold on;
axis equal;
% quiver(dmesh.r(1,:), dmesh.r(2,:), dmesh.n(1,:), dmesh.n(2,:),'k');
quiver(dmesh.r(1,:), dmesh.r(2,:), dmesh.pseudo_normals(1,:), dmesh.pseudo_normals(2,:),'r');