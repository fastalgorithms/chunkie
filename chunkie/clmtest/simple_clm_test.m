clear

% Setup geometry;
geom_class = [];
b = 4;
geom_class.xylim = [-b b -b b];
ndomain = 2;
geom_class.ndomain = ndomain;
ncurve = 1;
geom_class.ncurve = ncurve;

rn = ones(ndomain,1);
geom_class.rn = rn;
lambda = 0.01;
geom_class.lambda = lambda;
geom_class.mode = 'te';

a = 3;
verts = zeros(2,2);
verts(1,1) = -a;
verts(1,2) = a;

geom_class.verts = verts;

curves = cell(1,ncurve);
curves{1}.curvetype = 1;
curves{1}.vert_list = [1 2];
curves{1}.curve_id = 1;

geom_class.curves = curves;

clist = cell(1,ndomain);
clist{1} = [1];
clist{2} = [-1];

regions = cell(1,ndomain);
src = zeros(2,ndomain);
src(1,:) = [0.6, -0.9];
src(2,:) = [-3.1, 2.8];
for i=1:ndomain
    regions{i}.region_id = i;
    regions{i}.icurve_list = clist{i};
    regions{i}.is_inf = 0;
    if(i == 1)
        regions{i}.is_inf = 1;
    elseif (i==2)
        regions{i}.is_inf = -1;
    end
    regions{i}.src_in = src(:,i);
end
geom_class.regions = regions;
geom_class.lvert = 1;
geom_class.rvert = 2;
geom_class.ppw = 10;
clmparams = clm.setup_geom(geom_class);

[chnkr] = clm.get_geom_clmparams(clmparams);
figure(1)
clf
plot(chnkr,'k.')

chnkrtotal = merge(chnkr);


srcinfo = [];
targinfo = [];
srcinfo.r = clmparams.src_in(:,1);
targinfo.r = chnkrtotal.r(:,:);
rhs = chnk.helm2d.kern(clmparams.k(1),srcinfo,targinfo,'s');

fkern = @(s,t) chnk.helm2d.kern(clmparams.k(1),s,t,'d');
kvals = chunkerkernevalmat(chnkrtotal,fkern,clmparams.src_in(:,2));
pot = -2*kvals*rhs;

srcinfo = [];
srcinfo.r = clmparams.src_in(:,1);
targinfo = [];
targinfo.r = clmparams.src_in(:,2);
potex = chnk.helm2d.kern(clmparams.k(1),srcinfo,targinfo,'s');

fprintf('error in double layer test =%d\n',norm(pot-potex));



