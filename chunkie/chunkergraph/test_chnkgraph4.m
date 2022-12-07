%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple pentagonal chunkergraph 
%            and tests the Helmholtz transmission problem
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

zk_in = 2.0;
zk_ou = 1.0;

fkernd = @(s,t) chnk.helm2d.kern(zk_ou,s,t,'d')-...
    chnk.helm2d.kern(zk_in,s,t,'d');
fkerns = @(s,t) chnk.helm2d.kern(zk_ou,s,t,'s')-...
    chnk.helm2d.kern(zk_in,s,t,'s');
fkerndp = @(s,t) chnk.helm2d.kern(zk_ou,s,t,'dprime')-...
    chnk.helm2d.kern(zk_in,s,t,'dprime');
fkernsp = @(s,t) chnk.helm2d.kern(zk_ou,s,t,'sprime')-...
    chnk.helm2d.kern(zk_in,s,t,'sprime');

fkernmat = @(s,t,i,j)[fkernd(s,t), fkerns(s,t);fkerndp(s,t), fkernsp(s,t)];

opts = [];
[sysmat] = chunkgraphmat(cgrph.echnks,fkernmat,opts)

opts = [];
[sysmat] = chunkermat(cgrph,fkernmat,opts);
sysmat = sysmat - eye(size(sysmat,2))/2;

% generate some targets...

xs = -1:0.01:1;
ys = -1:0.01:1;
[X,Y] = meshgrid(xs,ys);
targs = [X(:).';Y(:).'];

srcinfo = [];
srcinfo.sources = cgrph.r(:,:);
w = weights(cgrph);
n = normals(cgrph);

% a quick hack to find the interior points

srcinfo.dipstr = w(:).';
srcinfo.dipvec = n(:,:); 
eps = 1E-8;
pg  = 0;
pgt = 1;
[U] = lfmm2d(eps,srcinfo,pg,targs,pgt);
U = U.pottarg;
inds = find(abs(U-2*pi)<pi/10);

%%%%%%%%%%%%%%%%%%
% generate the right hand side

x0 = 1.3;
y0 = 0.9;

srcinfo = [];
srcinfo.sources = cgrph.r(:,:);
w = weights(cgrph);
n = normals(cgrph);

s = [];
s.r = [x0;y0];
t = [];
t.r = cgrph.r(:,:);
rhs = chnk.helm2d.kern(zk,s,t,'s');
dens = sysmat\rhs;

srcinfo.dipstr = (w(:).*dens(:)).';
srcinfo.dipvec = n(:,:); 
fints = hfmm2d(eps,zk,srcinfo,pg,targs(:,inds),pgt);

t.r = targs(:,inds);
true_sol = chnk.helm2d.kern(zk,s,t,'s');

usol = zeros(size(targs,2),1);
uerr = usol;
usol(inds) = fints.pottarg;
uerr(inds) = usol(inds)-true_sol;
usol = reshape(usol,size(X));
uerr = reshape(uerr,size(X));