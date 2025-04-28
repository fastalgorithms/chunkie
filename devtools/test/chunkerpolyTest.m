chunkerpolyTest0();


function chunkerpolyTest0()
%CHUNKERPOLYTEST


% pre-defined vertices for a barbell shape

verts = chnk.demo.barbell(2.0,2.0,1.0,1.0);

% for barbell(2,2,1,1) area and length are

barb_area = 9;
barb_length = 16;

% rounded corner version

cparams = [];
cparams.widths = 0.1*ones(size(verts,2),1);
cparams.eps = 1e-8;

nv = size(verts,2);
edgevals = rand(3,nv); % constant values on edges to smooth out in arclength

p.k = 16; p.dim = 2;
chnkr = chunkerpoly(verts,cparams,p,edgevals);
chnkr = chnkr.sort();
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

figure(3)
clf
nchplot = 1:chnkr.nch;
x = chnkr.r(1,:,nchplot);
y =chnkr.r(2,:,nchplot);
z = chnkr.data(1,:,nchplot);
plot3(x(:),y(:),z(:))

% adaptive refinement in corners (no smoothing of edge data)

p.k = 16; p.dim = 2;
cparams = [];
cparams.rounded = false;
cparams.depth = 8;
chnkr2 = chunkerpoly(verts,cparams,p,edgevals);
chnkr2 = chnkr2.sort();
assert(checkadjinfo(chnkr2) == 0);

x = chnkr2.r(1,:); y = chnkr2.r(2,:); z = chnkr2.data(1,:);

figure(4)
clf
plot(chnkr2,'b-x')
hold on
quiver(chnkr2)

figure(5)
clf
plot3(x(:),y(:),z(:))

barb_area_2 = area(chnkr2);
err_area = abs(barb_area-barb_area_2)/abs(barb_area);
barb_length_2 = sum(sum(chnkr2.wts));
err_length = abs(barb_length -barb_length_2)/abs(barb_length);
fprintf('%5.2e : diff between true/computed area\n',err_area);
fprintf('%5.2e : diff between true/computed length\n',err_length);

%

verts = randn(2,5);

cparams = [];
cparams.widths = 0.1*ones(size(verts,2),1);
cparams.autowidths = true;
cparams.autowidthsfac = 0.1;
cparams.ifclosed = 0;
cparams.eps = 1e-3;

p.k = 16; p.dim = 2;
chnkr3 = chunkerpoly(verts,cparams,p);


figure(6)
clf
plot(chnkr3,'-x')
hold on
quiver(chnkr3)
axis equal

figure(7)
chnkr_ref3 = refine(chnkr3);
clf
plot(chnkr_ref3,'-x')
hold on
quiver(chnkr_ref3)
axis equal


end


