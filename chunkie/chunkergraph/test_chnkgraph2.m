%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple pentagonal chunkergraph 
%            and tests the interior Dirichlet Laplace problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

verts = exp(1i*2*pi*(0:4)/5);
verts = [real(verts);imag(verts)];

edge2verts = [-1, 1, 0, 0, 0; ...
               0,-1, 1, 0, 0; ...
               0, 0,-1, 1, 0; ...
               0, 0, 0,-1, 1; ...
               1, 0, 0, 0,-1];
edge2verts = sparse(edge2verts);


fchnks    = {};

prefs      = [];
prefs.chsmall = 1d-3;
[cgrph] = chunkgraphinit(verts,edge2verts,fchnks,prefs);

vstruc = procverts(cgrph);
rgns = findregions(cgrph);
cgrph = balance(cgrph);

fkern = @(s,t) chnk.lap2d.kern(s,t,'d');

opts = [];
[sysmat] = chunkermat(cgrph,fkern,opts);
sysmat = sysmat - eye(size(sysmat,2))/2;

rhs  = ones(size(cgrph.r(:,:),2),1);
dens = sysmat\rhs; 

% generate some targets...

xs = -1:0.01:1;
ys = -1:0.01:1;
[X,Y] = meshgrid(xs,ys);
targs = [X(:).';Y(:).'];

srcinfo = [];
srcinfo.sources = cgrph.r(:,:);
w = weights(cgrph);
n = normals(cgrph);

srcinfo.dipstr = w(:).';
srcinfo.dipvec = n(:,:); 
eps = 1E-8;
pg  = 0;
pgt = 1;
[U] = lfmm2d(eps,srcinfo,pg,targs,pgt);
U = U.pottarg;
inds = find(abs(U-2*pi)<pi/10);

srcinfo.dipstr = (w(:).*dens(:)).';
fints = lfmm2d(eps,srcinfo,pg,targs(:,inds),pgt);
usol = zeros(size(targs,2),1);
uerr = usol;
%%% differences in kernel convention
usol(inds) = -fints.pottarg/(2*pi);
uerr(inds) = abs(usol(inds)-1);
usol = reshape(usol,size(X));
uerr = reshape(uerr,size(X));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  now a harder problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x0 = 1.3;
y0 = 0.9;

s = [];
s.r = [x0;y0];
t = [];
t.r = cgrph.r(:,:);
rhs = chnk.lap2d.kern(s,t,'s');
dens2 = sysmat\rhs;

srcinfo.dipstr = (w(:).*dens2(:)).';
fints = lfmm2d(eps,srcinfo,pg,targs(:,inds),pgt);

t.r = targs(:,inds);
true_sol = chnk.lap2d.kern(s,t,'s');

usol2 = zeros(size(targs,2),1);
uerr2 = usol2;
%%% differences in kernel convention
usol2(inds) = -fints.pottarg/(2*pi);
uerr2(inds) = abs(usol2(inds)-true_sol);
usol2 = reshape(usol2,size(X));
uerr2 = reshape(uerr2,size(X));