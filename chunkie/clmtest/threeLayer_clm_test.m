clear;

% k1           s1 - s3
% ----- y = d
% k2    y = 0  s2 - s1
% ----- y = -d
% k1           s3 - s2

d = 1; dshift = d/2;
rn = [1.2,1,1.2];
src(1,:) = [0.6, 0.2,-0.9];
src(2,:) = [d+dshift,dshift,-d-dshift];


src_test = [0.3;dshift/2];



% get densities on each boundary
[clmparams1,chnkrtotal1] = setup_layer(rn(1:2),[src(:,2),src(:,1)],d);
[clmparams2,chnkrtotal2] = setup_layer(rn(2:3),[src(:,3),src(:,2)],-d);
k1 = clmparams1.k(1); k2 = clmparams1.k(2);

fkern_s_diff = @(s,t) -chnk.helm2d.kern(k1,s,t,'s') + ...
    chnk.helm2d.kern(k2,s,t,'s');
fkern_dprime_diff = @(s,t) chnk.helm2d.kern(k1,s,t,'dprime') - ...
    chnk.helm2d.kern(k2,s,t,'dprime');


fkern_s = @(s,t) -chnk.helm2d.kern(k2,s,t,'s');
fkern_d = @(s,t) chnk.helm2d.kern(k2,s,t,'d');
fkern_sprime = @(s,t) -chnk.helm2d.kern(k2,s,t,'sprime');
fkern_dprime = @(s,t) chnk.helm2d.kern(k2,s,t,'dprime');

opts = [];
opts.forcesmooth = true;


A11_21 = chunkerkernevalmat(chnkrtotal2,fkern_d,chnkrtotal1,opts);
A12_21 = chunkerkernevalmat(chnkrtotal2,fkern_s,chnkrtotal1,opts);
A21_21 = chunkerkernevalmat(chnkrtotal2,fkern_dprime,chnkrtotal1,opts);
A22_21 = chunkerkernevalmat(chnkrtotal2,fkern_sprime,chnkrtotal1,opts);



A11_12 = chunkerkernevalmat(chnkrtotal1,fkern_d,chnkrtotal2,opts);
A12_12 = chunkerkernevalmat(chnkrtotal1,fkern_s,chnkrtotal2,opts);
A21_12 = chunkerkernevalmat(chnkrtotal1,fkern_dprime,chnkrtotal2,opts);
A22_12 = chunkerkernevalmat(chnkrtotal1,fkern_sprime,chnkrtotal2,opts);



n = chnkrtotal1.npt;
tic, A12_1 = chunkermat(chnkrtotal1,fkern_s_diff); toc;
tic, A21_1 = chunkermat(chnkrtotal1,fkern_dprime_diff); toc;
A11 = -eye(n);
A22 = -eye(n);
Amat1 = [A11,A12_1;A21_1,A22];

A12_2 = -A12_1;
A21_2 = -A21_1;
Amat2 = [A11,A12_2;A21_2,A22];

Amat_total = [A11, A12_1, -A11_21, -A12_21; ...
              A21_1, A22, -A21_21, -A22_21;
              A11_12, A12_12, A11, A12_2;
              A21_12, A22_12, A21_2, A22];

% setup boundary data
rhs_1 = zeros(2*n,1); rhs_2 = zeros(2*n,1);

% (upper boundary)
[u1_1,grad1_1] = chnk.helm2d.green(k1,src(:,2),chnkrtotal1.r(:,:));
[u2_1,grad2_1] = chnk.helm2d.green(k2,src(:,3),chnkrtotal1.r(:,:));
% (lower boundary)
[u1_2,grad1_2] = chnk.helm2d.green(k2,src(:,3),chnkrtotal2.r(:,:));
[u2_2,grad2_2] = chnk.helm2d.green(k1,src(:,1),chnkrtotal2.r(:,:));


un1_1 = -grad1_1(:,:,2);
un2_1 = -grad2_1(:,:,2);

un1_2 = -grad1_2(:,:,2);
un2_2 = -grad2_2(:,:,2);


rhs_1(1:n) = (u1_1-u2_1);
rhs_1((n+1):end) = (un1_1-un2_1);

rhs_2(1:n) = (u1_2-u2_2);
rhs_2((n+1):end) = (un1_2-un2_2);


rhs_total = zeros(4*n,1);
rhs_total(1:2*n) = rhs_1;
rhs_total((2*n+1):end) = rhs_2;


ifgreenfun = 1;
if( ifgreenfun)
    [u2_1,grad2_1] = chnk.helm2d.green(k2,src(:,2),chnkrtotal1.r(:,:));
    [u1_2,grad1_2] = chnk.helm2d.green(k2,src(:,2),chnkrtotal2.r(:,:));
    rhs_total = zeros(4*n,1);
    rhs_total(1:n) = u2_1;
    rhs_total((n+1):(2*n)) = -grad2_1(:,:,2);
    
    rhs_total((2*n+1):(3*n)) = -u1_2;
    rhs_total((3*n+1):end) = grad1_2(:,:,2);
end

% solve (run gmres on dense(Amat) instead of backslash)

start = tic; b = gmres(Amat_total,rhs_total,[],1e-14,200); t1 = toc(start);

psi1 = b(1:n);
phi1 = b((n+1):2*n);

psi2 = b((2*n+1):3*n);
phi2 = b((3*n+1):end);

% Now test solution


fd1  = @(s,t) chnk.helm2d.kern(k1,s,t,'d');
fs1  = @(s,t) chnk.helm2d.kern(k1,s,t,'s');
fd2  = @(s,t) chnk.helm2d.kern(k2,s,t,'d');
fs2  = @(s,t) chnk.helm2d.kern(k2,s,t,'s');

psi_res1 = reshape(psi1,[chnkrtotal1.k,chnkrtotal1.nch]);
phi_res1 = reshape(phi1,[chnkrtotal1.k,chnkrtotal1.nch]);
psi_res2 = reshape(psi2,[chnkrtotal2.k,chnkrtotal2.nch]);
phi_res2 = reshape(phi2,[chnkrtotal2.k,chnkrtotal2.nch]);

uscat_d1 = chunkerkerneval(chnkrtotal1,fd1,psi_res1,src(:,1));
uscat_s1 = chunkerkerneval(chnkrtotal1,fs1,phi_res1,src(:,1));

uscat_d2 = chunkerkerneval(chnkrtotal1,fd2,psi_res1,src(:,2)) + ...
           chunkerkerneval(chnkrtotal2,fd2,psi_res2,src(:,2));
uscat_s2 = chunkerkerneval(chnkrtotal1,fs2,phi_res1,src(:,2)) + ...
           chunkerkerneval(chnkrtotal2,fs2,phi_res2,src(:,2));

uscat_d3 = chunkerkerneval(chnkrtotal2,fd1,psi_res2,src(:,3));
uscat_s3 = chunkerkerneval(chnkrtotal2,fs1,phi_res2,src(:,3));


sol1 = uscat_d1 - uscat_s1;
solex1 = chnk.helm2d.green(k1,src(:,1),src(:,2));
norm(sol1-solex1)

sol2 = uscat_d2 - uscat_s2;
solex2 = chnk.helm2d.green(k2,src(:,2),src(:,3));
norm(sol2-solex2)

sol3 = uscat_d3 - uscat_s3;
solex3 = chnk.helm2d.green(k1,src(:,3),src(:,1));
norm(sol3-solex3)


% at test point in region 2;
uscat_d2 = chunkerkerneval(chnkrtotal1,fd2,psi_res1,src_test) + ...
           chunkerkerneval(chnkrtotal2,fd2,psi_res2,src_test);
uscat_s2 = chunkerkerneval(chnkrtotal1,fs2,phi_res1,src_test) + ...
           chunkerkerneval(chnkrtotal2,fs2,phi_res2,src_test);
       
sol_test = uscat_d2 - uscat_s2;

sol_test_add = chnk.helm2d.green(k2,src(:,2),src_test);
sol_test_sub = chnk.helm2d.green(k1,src(:,2),src_test);

sol_test = sol_test + sol_test_add - sol_test_sub;


function [clmparams,chnkrtotal] = setup_layer(rn,src,layer_y)

% Setup geometry;
geom_class = [];
b = 4;
geom_class.xylim = [-b b -b b];
ndomain = 2;
geom_class.ndomain = ndomain;
ncurve = 1;
geom_class.ncurve = ncurve;


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
geom_class.ppw = 5;
clmparams = clm.setup_geom(geom_class);

[chnkr] = clm.get_geom_clmparams(clmparams);
%clf
%plot(chnkr,'k.')

chnkrtotal = merge(chnkr);
targs = chnkrtotal.r(:,:)+[0;layer_y];
chnkrtotal.r = reshape(targs,[2,chnkrtotal.k,chnkrtotal.nch]);
end
