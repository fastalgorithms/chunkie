nv = 10;
z = exp(1j*2*pi*(1:nv-1)/nv);
verts = [real(z); imag(z)];
nchs = randi([1,10], 1, nv);

umesh = [];
umesh.verts = verts;
umesh.pseudo_normals = get_pseudo_normal(verts);

qmesh = [];
qmesh
dmesh = [];
