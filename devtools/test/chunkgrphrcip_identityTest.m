% silence the near-singular warnings that occur in hypersingular corners
ws = warning('off','MATLAB:nearlySingularMatrix');

chunkgrphrcip_identityTest0();
chunkgrphrcip_identityTest1();
chunkgrphrcip_identityTest2();

% restore the previous warning state
warning(ws);


function chunkgrphrcip_identityTest0()
% constant rcip_identity c on a pentagon, checked against c times the
% default assembly

verts = exp(1i*2*pi*(0:4)/5);
verts = [real(verts);imag(verts)];

edge2verts = [-1, 1, 0, 0, 0; ...
               0,-1, 1, 0, 0; ...
               0, 0,-1, 1, 0; ...
               0, 0, 0,-1, 1; ...
               1, 0, 0, 0,-1];
edge2verts = sparse(edge2verts);
amp = 0.5;
frq = 6;
fchnks    = {};
for icurve = 1:size(edge2verts,1)
    fchnks{icurve} = @(t) sinearc(t,amp,frq);
end

prefs      = [];
[cgrph] = chunkgraph(verts,edge2verts,fchnks,prefs);


zk = 1.0;
fkern = @(s,t) -2*chnk.helm2d.kern(zk,s,t,'d');

c = 2.3;

opts = [];
opts.nonsmoothonly = false;
opts.rcip = true;

% default system I + chunkermat
opts.rcip_identity = 1;
[sysmat1] = chunkermat(cgrph,fkern,opts);
sys1 = eye(size(sysmat1,2)) + sysmat1;

% scaled system c*I + c*chunkermat = c*(I + chunkermat)
opts.rcip_identity = c;
fkernc = @(s,t) c*fkern(s,t);
[sysmatc] = chunkermat(cgrph,fkernc,opts);
sysc = c*eye(size(sysmatc,2)) + sysmatc;

err = norm(sysc-c*sys1)/norm(c*sys1);
fprintf('constant rcip_identity matrix discrepancy %5.2e\n',err);
assert(err < 1e-12);

end


function chunkgrphrcip_identityTest1()
% Per-edge rcip_identity on a two edge lens with mixed Dirichlet and Neumann.
% The edges share both corners, so a different scaling meets on each side of
% each corner. Checked against the default assembly with scaled kernels.

% two arcs between (-1,0) and (1,0), forming a closed lens
verts = [-1 1; 0 0];
edgesendverts = [1 1; 2 2];
fchnks = cell(1,2);
fchnks{1} = @(t) lensarc(t,0.5);  % bows up
fchnks{2} = @(t) lensarc(t,-0.5); % bows down

cparams = [];
cparams.nover = 2;
[cgrph] = chunkgraph(verts,edgesendverts,fchnks,cparams);

zk = 30;

edir = 1; % Dirichlet edge
eneu = 2; % Neumann edge

% scale so the diagonal is identity (required for rcip)
kerns(length([edir,eneu]),length([edir,eneu])) = kernel();
kerns(edir,edir) = -2*kernel('helm', 'd', zk);
kerns(eneu,edir) = -2*kernel('helm', 'dp', zk);

kerns(edir,eneu) = 2*kernel('helm', 's', zk);
kerns(eneu,eneu) = 2*kernel('helm', 'sp', zk);

nedges = length(cgrph.echnks);

rid = 1 + 0.5*(1:nedges).'; % different strength per edge

opts = [];
opts.rcip_identity = rid;
[sysmatv] = chunkermat(cgrph,kerns,opts);

% diagonal identity from the per-edge scalings
dval = zeros(cgrph.npt,1);
for ie = 1:nedges
    inds = cgrph.edgeinds(ie);
    dval(inds) = rid(ie);
end
sysv = diag(dval) + sysmatv;

% scale rows of edge ie by 1/rid(ie), assemble, then rescale by rid(ie)
kernscal = kerns;
for ie = 1:nedges
    for je = 1:nedges
        kernscal(ie,je) = kernscal(ie,je).*(1/rid(ie));
    end
end
[sysmat1] = chunkermat(cgrph,kernscal);
sys1 = eye(size(sysmat1,2)) + sysmat1;
sys1 = diag(dval)*sys1;

err = norm(sysv-sys1)/norm(sys1);
fprintf('per-edge rcip_identity matrix discrepancy %5.2e\n',err);
assert(err < 1e-12);

end


function chunkgrphrcip_identityTest2()
% The same lens as test1 but with opdims > 1. Checks the matrix and that
% rcipsav postprocessing matches the default


% two arcs between (-1,0) and (1,0), forming a closed lens
verts = [-1 1; 0 0];
edgesendverts = [1 1; 2 2];
fchnks = cell(1,2);
fchnks{1} = @(t) lensarc(t,0.5);  % bows up
fchnks{2} = @(t) lensarc(t,-0.5); % bows down

cparams = [];
cparams.nover = 2;
[cgrph] = chunkgraph(verts,edgesendverts,fchnks,cparams);

zk = 30;
nedges = length(cgrph.echnks);

% opdims 2 blocks, interleaved from scalar helmholtz kernels
sk = kernel('helm','s',zk);
dk = kernel('helm','d',zk);
spk = kernel('helm','sprime',zk);
dpk = kernel('helm','dprime',zk);
kblock = kernel([sk dk; spk dpk]);

rid = 1 + 0.5*(1:nedges).'; % different strength per edge

% scaled system, save full recursion for postprocessing
opts = [];
opts.rcip_identity = rid;
opts.rcip_savedepth = Inf;
[sysmatv,~,rcipsavv] = chunkermat(cgrph,kblock,opts);

% diagonal identity, ndim rows per point on each edge
ndim = kblock.opdims(1);
dval = [];
for ie = 1:nedges
    nrowe = cgrph.echnks(ie).npt*ndim;
    dval = [dval; rid(ie)*ones(nrowe,1)];
end
sysv = diag(dval) + sysmatv;

% scale rows of edge ie by 1/rid(ie), assemble, then rescale by rid(ie)
kernscal(nedges,nedges) = kernel();
for ie = 1:nedges
    for je = 1:nedges
        kernscal(ie,je) = kblock.*(1/rid(ie));
    end
end
opts1 = [];
opts1.rcip_savedepth = Inf;
[sysmat1,~,rcipsav1] = chunkermat(cgrph,kernscal,opts1);
sysscal = eye(size(sysmat1,2)) + sysmat1;
sys1 = diag(dval)*sysscal;

err = norm(sysv-sys1)/norm(sys1);
fprintf('per-edge opdims rcip_identity matrix discrepancy %5.2e\n',err);
assert(err < 1e-12);

% density is invariant: (D+K)x = b <=> (I+Dinv K)x = Dinv b
rhs = randn(size(sysv,1),1);
solv = sysv\rhs;
sol1 = sysscal\(diag(1./dval)*rhs);
assert(norm(solv-sol1)/norm(sol1) < 1e-10);

% postprocessing reconstructs the same fine density from either rcipsav
ndepth = 20;
for ivert = 1:length(rcipsavv)
    starind = rcipsavv{ivert}.starind;
    interpv = cell2mat(chnk.rcip.rhohatInterp(solv(starind),rcipsavv{ivert},ndepth));
    interp1 = cell2mat(chnk.rcip.rhohatInterp(sol1(starind),rcipsav1{ivert},ndepth));
    fprintf('postprocessing density discrepancy %5.2e\n', ...
        norm(interpv-interp1)/norm(interp1));
    assert(norm(interpv-interp1)/norm(interp1) < 1e-10);
end

end


function [r,d,d2] = sinearc(t,amp,frq)
xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end

function [r,d,d2] = lensarc(t,amp)
% arc from (-1,0) at t=0 to (1,0) at t=1, bowing to height amp at t=1/2
xs = -1 + 2*t;
ys = amp*sin(pi*t);
xp = 2*ones(size(t));
yp = amp*pi*cos(pi*t);
xpp = zeros(size(t));
ypp = -amp*pi*pi*sin(pi*t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
