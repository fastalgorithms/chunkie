nv = 3;
z = exp(1j*2*pi*(1:nv)/nv);
verts = [real(z); imag(z)];
nchs = 3*ones(nv,1);
k = 16;

kquad = 32;



umesh = chnk.smoother.get_umesh(verts);
dmesh = chnk.smoother.get_mesh(umesh, nchs, k);
qmesh = chnk.smoother.get_mesh(umesh, nchs, kquad);

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
quiver(dmesh.r(1,:), dmesh.r(2,:), dmesh.pseudo_normals(1,:), dmesh.pseudo_normals(2,:),'r');


[~, nd] = size(dmesh.r);
h = zeros(nd,1);
n_newton = 3;

sig0 = sqrt(5)*max(umesh.lengths);
lam = 10;


ww = qmesh.wts(:).';
dpx = dmesh.pseudo_normals(1,:).';
dpy = dmesh.pseudo_normals(2,:).';

rnx = qmesh.n(1,:).';
rny = qmesh.n(2,:).';

for i=1:n_newton
    [h] = chnk.smoother.newt_step(h,umesh,dmesh,qmesh,sig0,lam);
    % rt = dmesh.r;
    % rt(1,:) = rt(1,:) + (h.*dpx).';
    % rt(2,:) = rt(2,:) + (h.*dpy).';
    % 
    % [sig, sig_grad] = chnk.smoother.get_sigs(umesh, rt, sig0, lam);
    % [val, grad, hess, hess_sig] = chnk.smoother.green(qmesh.r, rt, sig);
    % gx = grad(:,:,1); 
    % gy = grad(:,:,2);
    % 
    % phi = -(gx*(rnx.*qmesh.wts) + gy*(rny.*qmesh.wts)) - 0.5;
    % 
    % 
    % h11 = hess(:,:,1,1).*ww;
    % h12 = hess(:,:,1,2).*ww;
    % h21 = hess(:,:,2,1).*ww;
    % h22 = hess(:,:,2,2).*ww;
    % 
    % h1sig = hess_sig(:,:,1).*ww;
    % h2sig = hess_sig(:,:,2).*ww;
    % 
    % dx1 = h11*rnx + h12*rny;
    % dy1 = h21*rnx + h22*rny;
    % 
    % dx2 = h1sig*rnx + h2sig*rny;
    % dy2 = dx2;
    % 
    % dsigx = sig_grad(:,1);
    % dsigy = sig_grad(:,2);
    % 
    % dx2 = dx2.*dsigx;
    % dy2 = dy2.*dsigy;
    % 
    % 
    % dphidh = (dx1 + dx2).*dpx + (dy1 + dy2).*dpy;
    % h = h - phi./dphidh;

    figure(i)    
    rt = dmesh.r;
    rt(1,:) = rt(1,:) + (h.*dpx).';
    rt(2,:) = rt(2,:) + (h.*dpy).';
    plot(rt(1,:), rt(2,:), 'k.')
 
end