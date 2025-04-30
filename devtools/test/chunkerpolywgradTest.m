%CHUNKERPOLYWGRADTEST

clearvars; clc;
addpaths_loc();

% pre-defined vertices for a barbell shape

verts = chnk.demo.barbell(2.0, 2.0, 1.0, 1.0);
nv = size(verts,2);

% for barbell(2,2,1,1) area and length are

barb_area = 9;
barb_length = 16;


%% linear functional test
%

verts = verts + 0.1*randn(2,nv);

fff = @(x,y) [cos(x.*y)+2*y,-sin(x.*y).*y,-sin(x.*y).*x+2*ones(size(x))];

rng(1);

pert = randn(2, nv);
pert = zeros(2, nv);
pert(1,1) = 1;
nt = 4;

p = [];
p.k = 16; p.dim = 2;
cparams = [];
cparams.rounded = true;
cparams.depth = 2;
cparams.autowidths = true;
cparams.smoothwidths = true;

[chnkr, igrad] = chunkerpoly(verts, cparams, p);
refopts = [];
% chnkr = refine(chnkr, refopts);


[~,wleg] = lege.exps(chnkr.k);
wtsraw = wleg(:)*ones(1, chnkr.nch); wtsraw = wtsraw(:);

x = chnkr.r(1,:); x = x(:);
y = chnkr.r(2,:); y = y(:);
wts = weights(chnkr);
fval = fff(x,y);
fint = sum(fval(:,1).*wts(:));
grad = vertgrad(chnkr, igrad, nv); 
g = grad*pert(:);
grads = vertdsdtgrad(chnkr, igrad, nv);
gs = grads*pert(:);
gx = g(1:2:end); gy = g(2:2:end);
gint = sum((fval(:,2).*gx(:)+fval(:,3).*gy(:)).*wts(:)) + ... 
    sum(fval(:,1).*gs.*wts(:));



for i = 1:nt
    eps = 10^(-(i));
    
    verts2 = verts + pert*eps;
    [chnkr2, igrad2] = chunkerpoly(verts2, cparams, p);
    chnkr2 = refine(chnkr2,refopts);
    wtsraw2 = wleg(:)*ones(1, chnkr2.nch); wtsraw2 = wtsraw2(:);

    x2 = chnkr2.r(1,:); x2 = x2(:);
    y2 = chnkr2.r(2,:); y2 = y2(:);
    wts2 = weights(chnkr2);
    fval2 = fff(x2,y2);
    fint2 = sum(fval2(:,1).*wts2(:));    
    grad2 = vertgrad(chnkr2, igrad2, nv);
    g2 = grad2*pert(:);
    grads2 = vertdsdtgrad(chnkr2, igrad2, nv);
    gs2 = grads2*pert(:);

    gx2 = g2(1:2:end); gy2 = g2(2:2:end);
    gint2 = sum((fval2(:,2).*gx2(:)+fval2(:,3).*gy2(:)).*wts2(:)) + ... 
        sum(fval2(:,1).*gs2.*wts2(:));
    err1 = abs(1-eps*0.5*(gint+gint2)/(fint2-fint));
    
    
    fprintf('%5.2e grad approx error with eps = %5.2e\n',err1,eps);
end

%% linear functional test (rounded)  
%


fff = @(x,y) [cos(x)+2*y,-sin(x),2*ones(size(x))];

pert = randn(2,nv);
nt = 8;

p = [];
p.k = 16; p.dim = 2;
cparams = [];
cparams.rounded = true;
cparams.widths= .1*ones(nv,1);
cparams.autowidths = true;
cparams.autowidthsfac = 0.1;
cparams.smoothwidths = true;

[chnkr, igrad] = chunkerpoly(verts, cparams, p);
chnkr.npt
refopts =[];
chnkr = refine(chnkr, refopts);
chnkr.npt

[x, wleg] = lege.exps(chnkr.k);
wtsraw = wleg(:)*ones(1,chnkr.nch); wtsraw = wtsraw(:);

x = chnkr.r(1,:); x = x(:);
y = chnkr.r(2,:); y = y(:);
wts = weights(chnkr);
fval = fff(x,y);
fint = sum(fval(:,1).*wts(:));
grad = vertgrad(chnkr, igrad, nv); 
g = grad*pert(:);
grads = vertdsdtgrad(chnkr, igrad, nv);
gs = grads*pert(:);
gx = g(1:2:end); gy = g(2:2:end);
gint = sum((fval(:,2).*gx(:)+fval(:,3).*gy(:)).*wts(:)) + ... 
    sum(fval(:,1).*gs.*wts(:));
rnrm = normals(chnkr);

for i = 1:nt
    eps = 10^(-(i));
    
    verts2 = verts + pert*eps;
    [chnkr2, igrad2] = chunkerpoly(verts2, cparams, p);
    chnkr2 = refine(chnkr2,refopts);
    
    wtsraw2 = wleg(:)*ones(1, chnkr2.nch); wtsraw2=wtsraw2(:);

    x2 = chnkr2.r(1,:); x2 = x2(:);
    y2 = chnkr2.r(2,:); y2 = y2(:);
    wts2 = weights(chnkr2);
    fval2 = fff(x2,y2);
    fint2 = sum(fval2(:,1).*wts2(:));    
    grad2 = vertgrad(chnkr2, igrad2, nv);
    g2 = grad2*pert(:);
    grads2 = vertdsdtgrad(chnkr2, igrad2, nv);
    gs2 = grads2*pert(:);

    gx2 = g2(1:2:end); gy2 = g2(2:2:end);    
    gint2 = sum((fval2(:,2).*gx2(:)+fval2(:,3).*gy2(:)).*wts2(:)) + ... 
        sum(fval2(:,1).*gs2.*wts2(:));
    err1 = abs(1-eps*0.5*(gint+gint2)/(fint2-fint));
    
    
    fprintf('%5.2e grad approx error with eps = %5.2e\n',err1,eps);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UNIT TESTS
%

%%
% Test parameter derivatives of rounded coordinates 

nt = 3;
hbell = 1/8;
a= 0.1;
r1 = verts(:,3); r2 = verts(:,4); r3 = verts(:,5);
[w2,dw2dr1,dw2dr2,dw2dr3] = smoothminwidth(r1,r2,r3,a);
ttest = linspace(-1,1,30); ttest=ttest(:);

%dw2dr1=zeros(2,1);dw2dr2=zeros(2,1);dw2dr3=zeros(2,1);

[r,d,~] = froundnew(ttest,r1,r2,r3,w2,hbell);
[rd1,rd2,rd3,sd1,sd2,sd3] = ...
    froundnewgrad(ttest,r1,r2,r3,w2,dw2dr1,dw2dr2,dw2dr3,hbell);

s = sqrt(sum(d.^2,1));

pert1 = randn(2,1); pert2 = randn(2,1); pert3 = randn(2,1);
g = reshape(pert1.'*rd1(:,:) + pert2.'*rd2(:,:) + pert3.'*rd3(:,:),2,[]);
sg = pert1.'*sd1 + pert2.'*sd2 + pert3.'*sd3;

for i= 1:nt
    eps = 10^(-i);
    r1i = r1+eps*pert1;
    r2i = r2+eps*pert2;
    r3i = r3+eps*pert3;
    [w2i,dw2dr1i,dw2dr2i,dw2dr3i] = smoothminwidth(r1i,r2i,r3i,a);
    %w2i=w2; dw2dr1i=zeros(2,1);dw2dr2i=zeros(2,1);dw2dr3i=zeros(2,1);
    [ri,di,d2i] = froundnew(ttest,r1i,r2i,r3i,w2i,hbell);
    
    [rd1i,rd2i,rd3i,sd1i,sd2i,sd3i] = ...
        froundnewgrad(ttest,r1i,r2i,r3i,w2i,dw2dr1i,dw2dr2i,dw2dr3i,hbell);
    si = sqrt(sum(di.^2,1));
    gi = reshape(pert1.'*rd1i(:,:) + pert2.'*rd2i(:,:) + pert3.'*rd3i(:,:),2,[]);
    sgi = pert1.'*sd1i + pert2.'*sd2i + pert3.'*sd3i;
    err1 = abs(1-eps*0.5*(gi+g)./(ri-r));
    serr1 = abs(1-eps*0.5*(sgi+sg)./(si-s));
    
    fprintf('%5.2e grad approx error with eps    = %5.2e\n',max(err1(:)),eps);
    fprintf('%5.2e ds grad approx error with eps = %5.2e\n',max(serr1(:)),eps);
end

%%
% test t derivatives of rounded coordinates

nt = 3;
hbell = 1/8;
a= 0.1;
r1 = verts(:,3)+0.01*randn(2,1); r2 = verts(:,4)+0.01*randn(2,1); r3 = verts(:,5)+0.01*randn(2,1);
[w2] = smoothminwidth(r1,r2,r3,a);
t0 = -1+2*rand();

[r,d,d2] = froundnew(t0,r1,r2,r3,w2,hbell);

pert = sign(randn());
g = pert*d;
g2 = pert*d2;

for i= 1:nt
    eps = 10^(-i);
    ti = t0+pert*eps;
    [ri,di,d2i] = froundnew(ti,r1,r2,r3,w2,hbell);
    gi = pert*di;
    g2i = pert*d2i;
    err1 = abs(1-eps*0.5*(gi+g)./(ri-r));
    err2 = abs(1-eps*0.5*(g2i+g2)./(di-d));
    fprintf('%5.2e der approx error with eps  = %5.2e\n',max(err1(:)),eps);
    fprintf('%5.2e der2 approx error with eps = %5.2e\n',max(err2(:)),eps);
end

%%
% test derivatives of the smoothmin width function

nt = 3;
hbell = 1/8;
a= 0.1;
r1 = verts(:,3); r2 = verts(:,4); r3 = verts(:,5);
[w2,dw2dr1,dw2dr2,dw2dr3] = smoothminwidth(r1,r2,r3,a);

pert1 = randn(2,1); pert2 = randn(2,1); pert3 = randn(2,1);
g = pert1.'*dw2dr1 + pert2.'*dw2dr2 + pert3.'*dw2dr3;

for i= 1:nt
    eps = 10^(-i);
    r1i = r1+eps*pert1;
    r2i = r2+eps*pert2;
    r3i = r3+eps*pert3;
    [w2i,dw2dr1i,dw2dr2i,dw2dr3i] = smoothminwidth(r1i,r2i,r3i,a);
    gi = pert1.'*dw2dr1i + pert2.'*dw2dr2i + pert3.'*dw2dr3i;
    err1 = abs(1-eps*0.5*(gi+g)/(w2i-w2));
    
    fprintf('%5.2e grad approx error with eps = %5.2e\n',err1,eps);
end

%%
% test derivatives of smoothmin

a = abs(randn()); b=abs(randn());
pert = abs(randn(2,1));

nt=6;

[f,ga,gb] = smoothmin(a,b);
g = [ga, gb]*pert;

for i = 1:nt
    eps = 10^(-i);
    a2 = pert(1)*eps+a;
    b2 = pert(2)*eps+b;
    
    [f2,ga2,gb2]=smoothmin(a2,b2);
    g2 = [ga2,gb2] * pert;
    err1 = abs(1-eps*0.5*(g+g2)/(f2-f));
    
    fprintf('%5.2e grad approx error with eps = %5.2e\n',err1,eps);
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Functions to debug below
% 

function [sm,dsmda,dsmdb] = smoothmin(a,b)
%SMOOTHMIN smooth minimum of two positive numbers
%
% smoothmin(a,b) = 1/log(exp(1/a)+exp(1/b)-1)
%
% has the properties:
%
% min(a,b)/2 < smoothmin(a,b) < min(a,b)
% smoothmin(a,b) < smoothmin(c,b) if a<c
% smoothmin(a,b) = smoothmin(b,a)

% 1/log(exp(1/a)+exp(1/b)-1) is the actual formula
% we avoid overflow below


abmin = min(a,b);
e1 = exp(-1./abmin);
ea = exp(1./a-1./abmin);
eb = exp(1./b-1./abmin);

logsafe = 1./abmin+log(ea+eb-e1);

sm = 1./logsafe;

% the derivative with respect to a is 
% 
% e^(1/a)/( log^2(e^(1/a)+e^(1/b)-1) * a^2 * (e^(1/a)+e^(1/b)-1) )

dsmda = 1.0./(logsafe.^2.*a.^2.*(1+exp(1./b-1./a)-exp(-1./a)));
dsmdb = 1.0./(logsafe.^2.*b.^2.*(1+exp(1./a-1./b)-exp(-1./b)));

end

function [w,dwdr1,dwdr2,dwdr3] = smoothminwidth(r1,r2,r3,a)

ll = sqrt(sum( (r1-r2).^2));
lr = sqrt(sum( (r3-r2).^2));

[sm,dsmdl,dsmdr] = smoothmin(ll,lr);

w = a*sm;
dwdr1 = a*dsmdl*(r1-r2)/ll;
dwdr3 = a*dsmdr*(r3-r2)/lr;
dwdr2 = a*(dsmdl*(r2-r1)/ll + dsmdr*(r2-r3)/lr);

end

function [r,d,d2] = froundnew(t,r1,r2,r3,w,h)
%
% Make the rounded corner of width w corresponding to
% the vertices r1, r2, r3 (with r2 in the middle)
%
% Input:
%
% t - parameter in [-1,1]
% r1 - first vertex
% r2 - middle vertex
% r3 - last vertex
% w - width of corner (first point should be 
%            r2 + w*(r1-r2)/|r1-r2|)
% h - std dev of Gaussian used in rounding
%      recommend h = 1/8 for ~16,15,14 digits 
%      in position,first,second derivative
%
% Output:
%
% r,d,d2 - position, derivative, second derivative
%           of corner rounded according to specified
%           Gaussian std dev and width of corner
%

dim = size(r1,1);
m1=1; b0=0;
[y,dy,d2y] = chnk.spcl.absconvgauss(t,m1,b0,h);

r12 = r1-r2; u12 = r12/norm(r12);
r32 = r3-r2; u32 = r32/norm(r32);

rhs = w*[u12, u32];
xx = zeros(dim,2);
xx(1:2,1:2) = [-1 1; 1 1];
amat = rhs/xx;

r = zeros(dim,length(t));
d = zeros(dim,length(t));
d2 = zeros(dim,length(t));

r(1,:) = t;
r(2,:) = y;
d(1,:) = 1.0;
d(2,:) = dy;
d2(1,:) = 0.0;
d2(2,:) = d2y;

r = r2+amat*r;
d = amat*d;
d2 = amat*d2;

end

function [drr1,drr2,drr3,dsr1,dsr2,dsr3] = froundnewgrad(t,r1,r2,r3,w,...
                                                     dwdr1,dwdr2,dwdr3,h)
%
% Get the derivatives of the position and arclength density 
% with respect to the vertices for a rounded corner.
%
% Input:
%
% t - parameter in [-1,1]
% r1 - first vertex
% r2 - middle vertex
% r3 - last vertex
% w - width of corner (first point should be 
%            r2 + w*(r1-r2)/|r1-r2|)
% dwdr1 - gradient of w w.r.t. r1
% dwdr2 - gradient of w w.r.t. r2
% dwdr3 - gradient of w w.r.t. r3
% h - std dev of Gaussian used in rounding
%      recommend h = 1/8 for ~16,15,14 digits 
%      in position,first,second derivative
%
% Output:
%
% drr1 - gradient of position with respect to r1
% drr2 - gradient of position with respect to r2
% drr3 - gradient of position with respect to r3
% dsr1 - gradient of arclength density with respect to r1
% dsr2 - gradient of arclength density with respect to r2
% dsr3 - gradient of arclength density with respect to r3
%
%

dim = size(r1,1);
m1=1; b0=0;

[y,dy,~] = chnk.spcl.absconvgauss(t,m1,b0,h);

r12 = r1-r2; nr12 = norm(r12); u12 = r12/nr12;
r32 = r3-r2; nr32 = norm(r32); u32 = r32/nr32;

xx = zeros(dim,2);
xx(1:2,1:2) = [-1 1; 1 1];
rhs = [u12, u32];
amat = w*(rhs/xx);

drr1 = zeros(dim,dim,length(t));
drr2 = zeros(dim,dim,length(t));
drr3 = zeros(dim,dim,length(t));
dsr1 = zeros(dim,length(t));
dsr2 = zeros(dim,length(t));
dsr3 = zeros(dim,length(t));

r = zeros(dim,length(t));
d = zeros(dim,length(t));

for i = 1:dim
    
    onei = zeros(dim,1); onei(i)=1;
    
    du12dr1i = onei/nr12 - r12(i)*r12/(nr12^3);
    du32dr1i = zeros(dim,1);
    dwdr1i = dwdr1(i);
    
    du12dr3i = zeros(dim,1);
    du32dr3i = onei/nr32 - r32(i)*r32/(nr32^3);
    dwdr3i = dwdr3(i);
    
    du12dr2i = -du12dr1i;
    du32dr2i = -du32dr3i;
    dwdr2i = dwdr2(i);
    
    rhsdr1i = w*[du12dr1i, du32dr1i] + dwdr1i*[u12, u32];
    rhsdr2i = w*[du12dr2i, du32dr2i] + dwdr2i*[u12, u32];
    rhsdr3i = w*[du12dr3i, du32dr3i] + dwdr3i*[u12, u32];    
    
    amatdr1i = rhsdr1i/xx;
    amatdr2i = rhsdr2i/xx;
    amatdr3i = rhsdr3i/xx;
    
    r(1,:) = t;
    r(2,:) = y;
    
    d(1,:) = 1.0;
    d(2,:) = dy;
    
    ad = amat*d;
    s = sqrt(sum(ad.^2,1));
        
    
    drr1(i,:,:) = reshape(amatdr1i*r,size(drr1(i,:,:)));
    drr2(i,:,:) = reshape(onei + amatdr2i*r,size(drr2(i,:,:)));
    drr3(i,:,:) = reshape(amatdr3i*r,size(drr2(i,:,:)));
    
    dsr1(i,:) = (sum(ad.*(amatdr1i*d),1))./s;
    dsr2(i,:) = (sum(ad.*(amatdr2i*d),1))./s;
    dsr3(i,:) = (sum(ad.*(amatdr3i*d),1))./s;
end


end


function grad = vertdsdtgrad(chnkr, igall, nvert)
%VERTDSDTGRAD gradient of chunk arclength density with respect to each 
% vertex
%
% Syntax: 
%   grad = vertgrad(chnkr)
%
% Input:
%
% chnkr - chunk polygon object
%
% Output:
%
% grad - is a sparse 
%           (chnkr.npts) x (chnkr.dim \cdot chnkr.nvert) 
%        matrix which gives the gradient of the chunker coordinates with 
%        respect to the vertex locations
%


dim = chnkr.dim;
npt = chnkr.npt;

grad = sparse(npt,dim*nvert);

itmp1 = [0,(1:(nvert-1))];
itmp2 = 1:nvert; 

isgrad = [dim^2*nvert+1;dim^2*nvert] + dim * [itmp1;itmp2];
isg = igall(isgrad);

for i = 1:nvert
    jshift = (i-1)*chnkr.dim;
    j1 = 1:dim; j1 = j1(:); 
    
    ipts = 1:chnkr.npt; ipts = ipts(:);
    ii = repmat(ipts.',dim,1); ii = ii(:);
    n1 = numel(ipts);
    jj = repmat(j1,n1,1); jj = jshift + jj(:);
    isgi = isg(1,i):isg(2,i);
    vv = chnkr.data(isgi,:,:); vv = vv(:);
    
    grad = grad + sparse(ii,jj,vv,npt,dim*nvert);    
    
end

end


function [grad] = vertgrad(chnkr, igall, nvert)
%VERTGRAD gradient of chunk polygon nodes with respect to each vertex
%
% Syntax: 
%   grad = vertgrad(chnkr)
%
% Input:
%
% chnkr - chunk polygon object
%
% Output:
%
% grad - is a sparse 
%           (chnkr.dim \cdot chnkr.npts) x (chnkr.dim \cdot chnkr.nvert) 
%        matrix which gives the gradient of the chunker coordinates with 
%        respect to the vertex locations
%

dim = chnkr.dim;
npt = chnkr.npt;

grad = sparse(dim*npt,dim*nvert);

itmp1 = [0,(1:(nvert-1))];
itmp2 = 1:nvert; 

igrad = [1;0] +dim^2* [itmp1;itmp2];
isgrad = [dim^2*nvert+1;dim^2*nvert] + dim * [itmp1;itmp2];
ig = igall(igrad); ig = reshape(ig,[],nvert);
isg = igall(isgrad);

for i = 1:nvert
    jshift = (i-1)*chnkr.dim;
    j1 = 1:dim; j1 = j1(:); 
    
    ipts = 1:chnkr.npt; ipts = ipts(:);
    i1 = repmat(j1.',dim,1); i1 = i1(:);
    ii = i1 + (ipts.' - 1)*dim; ii = ii(:);
    n1 = numel(ipts);
    jj = repmat(j1,n1*dim,1); jj = jshift + jj(:);
    igi = ig(1,i):ig(2,i); igi=igi(:);
    vv = chnkr.data(igi,:,:); vv = vv(:);
    
    grad = grad + sparse(ii,jj,vv,dim*npt,dim*nvert);
    
end

end
