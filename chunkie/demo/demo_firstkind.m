%DEMO FIRST KIND INTEGRAL EQUATION ON A POLYGON
%
% Solve the interior Helmholtz Dirichlet problem on a polygon with the
% first-kind representation u = S[sigma]. 
%
%% Define geometry
nv = 5;
R = 1;
ang = 2*pi*(0:nv-1)/nv + pi/2;
verts = R*[cos(ang); sin(ang)];
edges = [1:nv; [2:nv, 1]];

cgrph = chunkgraph(verts, edges, []);

%% Setup kernel and build system matrix
zk = 1.0;

S = kernel('helm', 's', zk) + kernel.eye(0);

opts = []; opts.rcip = true;
A = chunkermat(cgrph, S, opts);

%% Solve with data from an exterior point source
src = [2.5*R; 0.7*R];
srcinfo = struct('r', src);
rhs = chnk.helm2d.kern(zk, srcinfo, cgrph, 's');

sig = A \ rhs;

%% Check accuracy inside
targ = [0.1; -0.05];
u_t = chunkerkerneval(cgrph, S, sig, targ);
uex = chnk.helm2d.kern(zk, srcinfo, struct('r',targ), 's');
fprintf('firstkind: interior field relerr = %5.2e\n', abs(u_t - uex)/abs(uex));

x1 = linspace(-1.3*R, 1.3*R, 300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];

in = chunkerinterior(cgrph, {x1,x1});

uu = chunkerkerneval(cgrph, S, sig, targs);
uref = chnk.helm2d.kern(zk, srcinfo, struct('r',targs), 's');
u = nan(size(xx(:)));
u(in) = real(uu(in) - uref(in));

figure(1); clf
plot(cgrph, 'LineWidth', 2);
hold on
pcolor(xx, yy, reshape(log10(abs(u)+1e-16), size(xx))); shading interp
axis equal tight
colorbar
title('log_{10} interior error')