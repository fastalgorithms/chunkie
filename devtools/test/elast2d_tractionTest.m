elast2d_tractionTest0();
elast2d_tractionTest1();


function elast2d_tractionTest0()
%ELAST2D_TRACTIONTEST
%
% analytic solution test for the interior traction (Neumann) problem of
% 2D linear elasticity on a square, discretized as a chunkgraph with
% RCIP corner corrections. The boundary data is the traction of a sum of
% Kelvinlets (point forces) placed outside the domain. The solution is
% represented with the "string" single layer of kernel.elast2d_string,
% u = S_str[sigma]. The string attached to each boundary point y points
% radially away from the center of the square, (y - c0)/|y - c0|, so the
% string directions are continuous across the corners and never cross
% the (convex) domain. The 'strac' kernel does not depend on the string
% direction, so this choice only affects the evaluation of u (by a
% constant translation in the interior, checked below). The traction of
% S_str has no principal value part, and its interior limit is
%
%       t = sigma + mu*K sigma,   K = kernel.elast2d_string('strac',...)
%
% (the 'strac' kernel is the traction divided by mu). This is already of
% the form I + K needed by RCIP, so we solve
%
%       (I + mu*K) sigma = t.
%
% The interior traction problem only determines u up to a rigid motion,
% which is removed by a least squares fit before comparing. The log10
% error is plotted on a grid inside the geometry.
%
% With RCIP, the density on the corner panels is the compressed one,
% which is only accurate for targets well separated from the corner on
% the scale of the coarse panels, so the error is measured, the rigid
% motion fitted, and the assert applied on targets more than 0.3 from
% the corners; the plot shows the error everywhere.

iseed = 8675309;
rng(iseed,'twister');

% Lame parameters

lam = 1.5;
mu = 2.1;

% geometry: square [-1,1]^2 as a chunkgraph (counterclockwise, so the
% normals point out of the domain)

verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];
cparams = []; cparams.nchmin = 8;
pref = []; pref.k = 16;
start = tic; cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);
t1 = toc(start);
fprintf('%5.2e s : time to build geo, nch = %d, npt = %d\n',t1, ...
    cgrph.nch,cgrph.npt);

% point forces outside the domain

ns = 6;
ts = 2*pi*rand(1,ns);
src = []; src.r = 2.5*[cos(ts); sin(ts)]; src.n = randn(2,ns);
f = randn(2*ns,1);

% Kelvinlet kernels (exact solution) and string kernels (representation)

kerns = kernel.elast2d('s',lam,mu);
kernstrac = kernel.elast2d('strac',lam,mu);

ifimage = false;
kernstrac_str = kernel.elast2d_string('strac',lam,mu,ifimage);

% strings along the outward normal (default) and bent radially away from
% the center c0

c0 = [0;0];
kerns_strn = kernel.elast2d_string('s',lam,mu,ifimage);
kerns_str = kerns_strn;
kerns_str.eval = @(s,t) chnk.elast2d_string.kern(lam,mu, ...
    radialstrings(s,c0),t,'s',ifimage);

% traction data on the boundary

bdry = []; bdry.r = cgrph.r(:,:); bdry.n = cgrph.n(:,:);
bdry.d = cgrph.d(:,:);
tracbdry = kernstrac.eval(src,bdry)*f;

% build system (with RCIP at the corners) and solve

opts = []; opts.rcip = true;
start = tic; sysmat = chunkermat(cgrph,mu*kernstrac_str,opts); 
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1);

sys = eye(2*cgrph.npt) + sysmat;
rhs = tracbdry;
start = tic; [sigma,~,relres,iter] = gmres(sys,rhs,[],1e-13,200);
t1 = toc(start);
fprintf('%5.2e s : time for gmres, %d iters, relres %5.2e\n',t1, ...
    iter(end),relres);

% targets on a grid inside the geometry (boundary excluded)

nplot = 200;
rmin = min(cgrph.r(:,:),[],2); rmax = max(cgrph.r(:,:),[],2);
xs = linspace(rmin(1),rmax(1),nplot+2); xs = xs(2:end-1);
ys = linspace(rmin(2),rmax(2),nplot+2); ys = ys(2:end-1);
[xx,yy] = meshgrid(xs,ys);
targs = [xx(:).'; yy(:).'];

in = chunkerinterior(cgrph,targs);
targin = targs(:,in);
nt = size(targin,2);

% exact and computed displacements

targinfo = []; targinfo.r = targin;
uex = kerns.eval(src,targinfo)*f;

start = tic; usol = chunkerkerneval(cgrph,kerns_str,sigma,targin); t1 = toc(start);
fprintf('%5.2e s : time to evaluate at %d targets\n',t1,nt);
usol = usol(:);

% bending the strings should only change u by a translation: compare with
% normal strings on a few targets

isub = 1:max(1,floor(nt/50)):nt;
usub = reshape(usol,2,nt); usub = usub(:,isub);
usubn = chunkerkerneval(cgrph,kerns_strn,sigma,targin(:,isub));
dsub = usub - reshape(usubn,2,[]);
errtrans = max(vecnorm(dsub - mean(dsub,2)))/max(vecnorm(usub));
fprintf('radial vs normal strings: translation %5.2e %5.2e, ', ...
    mean(dsub,2));
fprintf('non-translation part %5.2e\n',errtrans);

% targets away from the corners
dcorner = min(vecnorm(targin - permute(verts,[1 3 2]),2,1),[],3);
far = dcorner > 0.3;

% remove rigid motion (translations + rotation) by least squares
% (fitted away from the corners)

rotfld = [targin(2,:); -targin(1,:)];
rmat = [repmat(eye(2),nt,1), rotfld(:)];
ifar = reshape([far; far],[],1);
cfs = rmat(ifar,:)\(uex(ifar)-usol(ifar));
usol = usol + rmat*cfs;

% errors

err = reshape(uex-usol,2,nt);
unorm = max(vecnorm(reshape(uex,2,nt)));
errpt = vecnorm(err)/unorm;
relerr = norm(err(:))/norm(uex);
fprintf('all targets: relative l2 error %5.2e, max pointwise rel. error %5.2e\n', ...
    relerr,max(errpt));
errfar = max(errpt(far));
fprintf('targets > 0.3 from corners: max pointwise rel. error %5.2e\n', ...
    errfar);

% plot the computed solution (with the rigid motion removed) and the
% log10 error

usol2 = reshape(usol,2,nt);
fplot = {usol2(1,:), usol2(2,:), vecnorm(usol2), log10(errpt)};
tlabs = {'$u_1$', '$u_2$', '$|u|$', '$\log_{10}$(error)'};

figure(1); clf
tiledlayout(2,2,'TileSpacing','compact');
for ip = 1:4
    nexttile
    zz = nan(size(xx));
    zz(in) = fplot{ip};
    h = pcolor(xx,yy,zz);
    set(h,'EdgeColor','none')
    hold on
    plot(cgrph,'k','LineWidth',2)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
sgtitle('interior traction problem on a square (strings)', ...
    'Interpreter','latex')

assert(errfar < 1e-8);
assert(errtrans < 1e-12);

end



function elast2d_tractionTest1()
%ELAST2D_TRACTIONTEST1
%
% the same traction problem with the data of the demo: traction (1,1)
% on the top edge, 0 on the other three edges. There is no exact
% solution, so two checks are made, both away from the corners (where
% the data is discontinuous and the RCIP density is the compressed one):
%
% (a) refinement: the solutions with nchmin = 8 and nchmin = 32 agree on
%     a grid of interior targets more than 0.3 from the corners;
% (b) boundary condition (nchmin = 8): the traction mu*Kstrac[sigma] is evaluated at x - d*n,
%     d = h, 2h, ..., 5h (h = 1e-3), for points x on each edge with
%     edge parameter in [-0.7,0.7], extrapolated to d = 0 (error O(h^5)),
%     and compared with the data.
%
% The data (1,1) on top has net force (2,2) and torque -2, so its
% L2(boundary) projection onto the rigid motions is removed first; in
% closed form this is (1,1)*[top] - (0.25,0.25) + 0.1875*(-y,x), the same
% as in demo_elast2d_square_traction. The rigid motion of u is fixed by
% making u orthogonal to the rigid motions on the targets.

lam = 1.5;
mu = 2.1;
kernstrac_str = kernel.elast2d_string('strac',lam,mu,false);
kerns_str = kernel.elast2d_string('s',lam,mu,false);
verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];

dat = @(x) [1;1]*(abs(x(2,:)-1) < 1e-12) - 0.25 + ...
    0.1875*[-x(2,:); x(1,:)];

% interior targets more than 0.3 from the corners

xs = linspace(-1,1,52); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
dcorner = min(vecnorm(targs - permute(verts,[1 3 2]),2,1),[],3);
targs = targs(:,dcorner > 0.3);

% boundary check points and inward offsets

sb = linspace(-0.7,0.7,15); e1 = ones(size(sb)); e0 = zeros(size(sb));
xb = [sb, e1, sb, -e1; -e1, sb, e1, sb];
nb = [e0, e1, e0, -e1; -e1, e0, e1, e0];
fb = dat(xb);
h = 1e-3;
wex = [5 -10 10 -5 1];

nchs = [8 32];
u = cell(2,1);
for ich = 1:2
    cparams = []; cparams.nchmin = nchs(ich);
    pref = []; pref.k = 16;
    cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);
    rhs = dat(cgrph.r(:,:));
    opts = []; opts.rcip = true;
    sys = eye(2*cgrph.npt) + chunkermat(cgrph,mu*kernstrac_str,opts);
    [sigma,~,relres,iter] = gmres(sys,rhs(:),[],1e-13,200);
    fprintf('top (1,1) data, nchmin = %2d: gmres %d iters, relres %5.2e\n', ...
        nchs(ich),iter(end),relres);
    ui = reshape(chunkerkerneval(cgrph,kerns_str,sigma,targs),2,[]);
    % fix the rigid motion: no translation/rotation component on targs
    rm = [repmat(eye(2),size(targs,2),1), ...
        reshape([-targs(2,:); targs(1,:)],[],1)];
    ui = reshape(ui(:) - rm*(rm\ui(:)),2,[]);
    u{ich} = ui;
    if ich == 1
        cg1 = cgrph; sig1 = sigma;
    end
end

% (a) refinement

errref = max(vecnorm(u{1}-u{2}))/max(vecnorm(u{2}));
fprintf('top (1,1) data: nchmin 8 vs 32, max rel. difference %5.2e\n', ...
    errref);

% (b) boundary condition (on the nchmin = 8 discretization)

bc = zeros(size(xb));
for id = 1:5
    t = []; t.r = xb - id*h*nb; t.n = nb;
    v = mu*chunkerkerneval(cg1,kernstrac_str,sig1,t);
    bc = bc + wex(id)*reshape(v,2,[]);
end
errbc = max(vecnorm(bc - fb))/max(vecnorm(fb));
fprintf('top (1,1) data: boundary condition max rel. error %5.2e\n', ...
    errbc);

% plot the solution

nplot = 100;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
tp = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,tp);
up = reshape(chunkerkerneval(cgrph,kerns_str,sigma,tp(:,in)),2,[]);
rm = [repmat(eye(2),size(up,2),1), ...
    reshape([-tp(2,in); tp(1,in)],[],1)];
up = reshape(up(:) - rm*(rm\up(:)),2,[]);
fplot = {up(1,:), up(2,:), vecnorm(up)};
tlabs = {'$u_1$', '$u_2$', '$|u|$'};
figure(2); clf
for ip = 1:3
    subplot(1,3,ip)
    zz = nan(size(xx));
    zz(in) = fplot{ip};
    hh = pcolor(xx,yy,zz);
    set(hh,'EdgeColor','none')
    hold on
    plot(cgrph,'k','LineWidth',2)
    plot([-1 1],[1 1],'r','LineWidth',3)
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
if exist('sgtitle','file')
    sgtitle('traction $(1,1)$ on top (red), $0$ elsewhere, rigid parts removed','Interpreter','latex')
end

assert(errref < 2e-8);
assert(errbc < 1e-9);

end

function s = radialstrings(s,c0)
% replace the source normals (which set the string directions in
% chnk.elast2d_string.kern) by unit vectors pointing away from c0
dr = s.r(:,:) - c0(:);
s.n = reshape(dr./vecnorm(dr),size(s.r));
end
