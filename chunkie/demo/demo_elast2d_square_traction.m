%DEMO_ELAST2D_SQUARE_TRACTION
%
% Interior traction (Neumann) problem of 2D linear elasticity on the
% square [-1,1]^2, with traction (1,1) prescribed on the top edge and
% zero traction on the other three edges.
%
% As given, this data is not compatible: the interior traction problem
% requires zero net force and torque, but here the net force is (2,2)
% and the torque about the origin is -2. We make it compatible by
% subtracting the L2(boundary) projection of the data onto the rigid
% motions (two translations and one rotation). The other edges then
% carry a small uniform/rotational traction.
%
% The solution is represented with the "string" single layer of
% kernel.elast2d_string, u = S_str[sigma], with the strings pointing
% radially away from the center of the square, and the corners are
% handled with RCIP. The traction of S_str has interior limit
% sigma + mu*K sigma (K = 'strac' string kernel, which is traction/mu),
% so we solve
%
%       (I + mu*K) sigma = t.
%
% u is only determined up to a rigid motion; the plotted solution is
% normalized to have no rigid-motion component (least squares on the
% plotting grid).

% BEGIN DEMO ELAST2D SQUARE TRACTION
%% Lame parameters

lam = 1.5;
mu = 2.1;

%% geometry: square as a chunkgraph (counterclockwise, outward normals)

verts = [-1 1 1 -1; -1 -1 1 1];
edgesendverts = [1:4; [2:4 1]];   % edge 3 (vertex 3 -> 4) is the top
cparams = []; cparams.nchmin = 8;
pref = []; pref.k = 16;
cgrph = chunkgraph(verts,edgesendverts,[],cparams,pref);

%% boundary data: (1,1) on the top edge, zero elsewhere

rb = cgrph.r(:,:);
npt = cgrph.npt;
itop = abs(rb(2,:) - 1) < 1e-12;
trac = zeros(2,npt);
trac(:,itop) = 1;
trac = trac(:);

% project out the rigid motions in L2(boundary)

wts = cgrph.wts(:);
w2 = repmat(wts.',2,1); w2 = w2(:);
rot = [-rb(2,:); rb(1,:)];
rmat = [repmat(eye(2),npt,1), rot(:)];
cfs = (rmat.'*(w2.*rmat)) \ (rmat.'*(w2.*trac));
trac = trac - rmat*cfs;

fprintf('rigid part removed: force (%5.2e, %5.2e), rotation %5.2e\n',cfs);
fprintf('compatibility check: net force (%5.2e, %5.2e), torque %5.2e\n', ...
    rmat.'*(w2.*trac));

%% kernels

ifimage = false;
c0 = [0;0];
kernstrac_str = kernel.elast2d_string('strac',lam,mu,ifimage);
kerns_str = kernel.elast2d_string('s',lam,mu,ifimage);
kerns_str.eval = @(s,t) chnk.elast2d_string.kern(lam,mu, ...
    radialstrings(s,c0),t,'s',ifimage);

%% build system (with RCIP at the corners) and solve

opts = []; opts.rcip = true;
start = tic; sysmat = chunkermat(cgrph,mu*kernstrac_str,opts);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1);

sys = eye(2*npt) + sysmat;
start = tic; [sigma,~,relres,iter] = gmres(sys,trac,[],1e-13,200);
t1 = toc(start);
fprintf('%5.2e s : time for gmres, %d iters, relres %5.2e\n',t1, ...
    iter(end),relres);

%% evaluate on a grid inside the square

nplot = 200;
xs = linspace(-1,1,nplot+2); xs = xs(2:end-1);
[xx,yy] = meshgrid(xs,xs);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(cgrph,targs);
targin = targs(:,in);
nt = size(targin,2);

start = tic; u = chunkerkerneval(cgrph,kerns_str,sigma,targin);
t1 = toc(start);
fprintf('%5.2e s : time to evaluate at %d targets\n',t1,nt);
u = u(:);

% normalize: remove the rigid-motion component on the grid

rotg = [-targin(2,:); targin(1,:)];
rmatg = [repmat(eye(2),nt,1), rotg(:)];
u = u - rmatg*(rmatg\u);
u = reshape(u,2,nt);

%% plot

fplot = {u(1,:), u(2,:), vecnorm(u)};
tlabs = {'$u_1$', '$u_2$', '$|u|$'};

figure(1); clf
tiledlayout(1,3,'TileSpacing','compact');
for ip = 1:3
    nexttile
    zz = nan(size(xx));
    zz(in) = fplot{ip};
    h = pcolor(xx,yy,zz);
    set(h,'EdgeColor','none')
    hold on
    plot(cgrph,'k','LineWidth',2)
    if ip == 3
        u1 = nan(size(xx)); u1(in) = u(1,:);
        u2 = nan(size(xx)); u2(in) = u(2,:);
        isk = 1:12:nplot;
        quiver(xx(isk,isk),yy(isk,isk),u1(isk,isk),u2(isk,isk),'w');
    end
    axis equal tight
    set(gca,"box","off","Xtick",[],"Ytick",[]);
    colorbar()
    title(tlabs{ip},'Interpreter','latex','FontSize',12)
end
sgtitle(['traction $(1,1)$ on top, zero elsewhere ' ...
    '(rigid part projected out)'],'Interpreter','latex')
% END DEMO ELAST2D SQUARE TRACTION


function s = radialstrings(s,c0)
% replace the source normals (which set the string directions in
% chnk.elast2d_string.kern) by unit vectors pointing away from c0
dr = s.r(:,:) - c0(:);
s.n = reshape(dr./vecnorm(dr),size(s.r));
end
