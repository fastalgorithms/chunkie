%TEST_CHUNKPOLY

addpaths_loc();


verts = chnk.demo.barbell(2.0,2.0,1.0,1.0);

cparams = [];
cparams.widths = 0.1*ones(size(verts,2),1);
cparams.eps = 1e-3;

p.k = 16; p.dim = 2;
chnkr = chunkpoly(verts,cparams,p);

assert(checkadjinfo(chnkr) == 0);

figure(1)
clf
plot(chnkr,'-x')
hold on
quiver(chnkr)
axis equal

figure(2)
chnkr_ref = refine(chnkr);
clf
plot(chnkr_ref,'-x')
hold on
quiver(chnkr_ref)
axis equal

%

verts = randn(2,5);

cparams = [];
cparams.widths = 0.1*ones(size(verts,2),1);
cparams.autowidths = true;
cparams.autowidthsfac = 0.1;
cparams.ifclosed = 0;
cparams.eps = 1e-3;

p.k = 16; p.dim = 2;
chnkr2 = chunkpoly(verts,cparams,p);


figure(3)
clf
plot(chnkr2,'-x')
hold on
quiver(chnkr2)
axis equal

figure(4)
chnkr_ref2 = refine(chnkr2);
clf
plot(chnkr_ref2,'-x')
hold on
quiver(chnkr_ref2)
axis equal
