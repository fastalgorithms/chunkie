clear

% Setup geometry;
geom_class = [];
b = 4;
geom_class.xylim = [-b b -b b];
ndomain = 2;
geom_class.ndomain = ndomain;
ncurve = 1;
geom_class.ncurve = ncurve;

rn = [1,1.2];
geom_class.rn = rn;
lambda = 2.0;
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

%% Now test transmission problem
% setup matrices
fkern_s_diff = @(s,t) -chnk.helm2d.kern(clmparams.k(1),s,t,'s') + ...
    chnk.helm2d.kern(clmparams.k(2),s,t,'s');
fkern_dprime_diff = @(s,t) chnk.helm2d.kern(clmparams.k(1),s,t,'dprime') - ...
    chnk.helm2d.kern(clmparams.k(2),s,t,'dprime');

n = chnkrtotal.npt;
tic, A12 = chunkermat(chnkrtotal,fkern_s_diff); toc;
tic, A21 = chunkermat(chnkrtotal,fkern_dprime_diff); toc;
A11 = -eye(n);
A22 = -eye(n);
Amat = [A11,A12;A21,A22];

% setup boundary data
rhs = zeros(2*n,1);

[u1,grad1] = chnk.helm2d.green(clmparams.k(1),src(:,1),chnkrtotal.r(:,:));
[u2,grad2] = chnk.helm2d.green(clmparams.k(2),src(:,2),chnkrtotal.r(:,:));

un1 = -grad1(:,:,2);
un2 = -grad2(:,:,2);

rhs(1:n) = (u1-u2);
rhs((n+1):end) = (un1-un2);

% solve
b = Amat\rhs;
psi = b(1:n);
phi = b((n+1):end);

% Now test solution
k1 = clmparams.k(1);
k2 = clmparams.k(2);

fd1  = @(s,t) chnk.helm2d.kern(k1,s,t,'d');
fs1  = @(s,t) chnk.helm2d.kern(k1,s,t,'s');

fd2  = @(s,t) chnk.helm2d.kern(k2,s,t,'d');
fs2  = @(s,t) chnk.helm2d.kern(k2,s,t,'s');

psi_res = reshape(psi,[chnkrtotal.k,chnkrtotal.nch]);
phi_res = reshape(phi,[chnkrtotal.k,chnkrtotal.nch]);

uscat_d1 = chunkerkerneval(chnkrtotal,fd1,psi_res,src(:,2));
uscat_s1 = chunkerkerneval(chnkrtotal,fs1,phi_res,src(:,2));

sol1 = uscat_d1 - uscat_s1;
solex1 = chnk.helm2d.green(k1,src(:,1),src(:,2));
err1 = norm(sol1-solex1);
fprintf('error in potential in region1 = %d\n',err1);

uscat_d2 = chunkerkerneval(chnkrtotal,fd2,psi_res,src(:,1));
uscat_s2 = chunkerkerneval(chnkrtotal,fs2,phi_res,src(:,1));

sol2 = uscat_d2 - uscat_s2;
solex2 = chnk.helm2d.green(k2,src(:,1),src(:,2));
err2 = norm(sol2-solex2);
fprintf('error in potential in region2 = %d\n',err2);