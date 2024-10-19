nv = 3;
z = exp(1j*2*pi*(1:nv)/nv);
verts = [real(z); imag(z)];
nchs = 3*ones(nv,1);
k = 16;

kquad = 32;



umesh = chnk.smoother.get_umesh(verts);
dmesh = chnk.smoother.get_mesh(umesh, nchs, k);
qmesh = chnk.smoother.get_mesh(umesh, nchs, kquad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nv = 3;
z = exp(1j*2*pi*(1:nv)/nv)+(0.9+1*1i);
verts2 = [real(z); imag(z)];
nchs = 3*ones(nv,1);
k = 16;

kquad = 32;



umesh2 = chnk.smoother.get_umesh(verts2);
dmesh2 = chnk.smoother.get_mesh(umesh2, nchs, k);
qmesh2 = chnk.smoother.get_mesh(umesh2, nchs, kquad);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

umeshes = {umesh,umesh2};
dmeshes = {dmesh,dmesh2};
qmeshes = {qmesh,qmesh2};
src_codes = [1,-1];
targ_codes = [0.5,0.5]*0;
[umesh,dmesh,qmesh,scales,levels] = chnk.smoother.merge(umeshes,dmeshes,...
    qmeshes,src_codes,targ_codes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
clf
verts = umesh.verts;
plot(verts(1,:), verts(2,:), 'k.','MarkerSize',20); hold on;
plot(umesh.centroids(1,:), umesh.centroids(2,:), 'r.', 'MarkerSize',20);
quiver(verts(1,:), verts(2,:), umesh.pseudo_normals(1,:), umesh.pseudo_normals(2,:),'k')
quiver(umesh.centroids(1,:), umesh.centroids(2,:), umesh.face_normals(1,:), umesh.face_normals(2,:),'r')

figure(2)
clf
plot(dmesh.r(1,:), dmesh.r(2,:),'k.'); hold on;
axis equal;
quiver(dmesh.r(1,:), dmesh.r(2,:), dmesh.pseudo_normals(1,:), dmesh.pseudo_normals(2,:),'r');


[~, nd] = size(dmesh.r);
h = zeros(nd,1);
n_newton = 20;

sig0 = sqrt(5)*max(umesh.lengths);
lam = 10;


ww = qmesh.wts(:).';
dpx = dmesh.pseudo_normals(1,:).';
dpy = dmesh.pseudo_normals(2,:).';

rnx = qmesh.n(1,:).';
rny = qmesh.n(2,:).';
opts = [];
opts.scales = scales;
opts.levels  = levels;
opts.step_fact = 0.9;
for i=1:n_newton
    [h] = chnk.smoother.newt_step(h,umesh,dmesh,qmesh,sig0,lam,opts);
    figure(i)    
    rt = dmesh.r;
    rt(1,:) = rt(1,:) + (h.*dpx).';
    rt(2,:) = rt(2,:) + (h.*dpy).';
    plot(rt(1,:), rt(2,:), 'k.')
 
end