%DEMO_CONCENTRIC_DOMAIN
%
% demonstrate Helmholtz transmission helper
% demonstrate chunkgraph regions
%
% code takes about 2 minutes


% wavenumbers
zks = 2.^(1+(1:5));

%% Build geometry
% define verts
verts_out = [[1;1],[1;-1],[-1;-1],[-1;1]];
verts_in = [[0;1],[1;0],[0;-1],[-1;0]];
verts_in2 = [[1/2;0],[-1/2;0]];
verts_in3 = verts_in2/2;

verts = [verts_out, verts_in, verts_in2, verts_in3];

% geometry will be a series of loops through subsets of vertices
id_vert_out = [5,1,6,2,7,3,8,4];
% call helper function
e2v_out = build_loop_edges(id_vert_out);

id_vert_in = [5:8];
e2v_in = build_loop_edges(id_vert_in);

id_vert_in2 = [5,9,7,10];
e2v_in2 =  build_loop_edges(id_vert_in2);

id_vert_in3 = [5,11,7,12];
e2v_in3 =  build_loop_edges(id_vert_in3);

% combine lists of edges
edge_2_verts = [e2v_out,e2v_in,e2v_in2,e2v_in3];

cparams = [];
cparams.maxchunklen = 4.0/max(zks);
cgrph = chunkgraph(verts,edge_2_verts,[],cparams);
figure(1);clf
% plot(cgrph)
% quiver(cgrph)
plot_regions(cgrph)

%% Build system
nreg = length(cgrph.regions);

% assign region wavenumbers
ks = [zks(1),repmat(zks(2),1,4),repmat(zks(3),1,2),repmat(zks(4),1,2),zks(5)];

% jump condition coefficients
coefs = ones(1,nreg);

% identify region on each side of each edge
edge_regs = zeros(2,size(edge_2_verts,2));
for i = 1:nreg
    id_edge = cgrph.regions{i}{1};
    edge_regs(1,id_edge(id_edge>0)) = i;
    edge_regs(2,-id_edge(id_edge<0)) = i;
end

% set up parameters for planewave data
opts = [];
opts.bdry_data_type = 'pw';
opts.exposed_curves = (edge_regs(1,:)==1) -(edge_regs(2,:)==1);

% build system and get boundary data using the transmission helper
[kerns, bdry_data, kerns_eval] = chnk.helm2d.transmission_helper(cgrph, ...
                                   ks, edge_regs, coefs, opts);
% build system matrix
tic;
[sysmat] = chunkermat(cgrph, kerns);
sysmat = sysmat + eye(size(sysmat,2));
tbuild = toc

% solve
tic; 
dens = sysmat\bdry_data;
tsolve = toc

%% Compute field

L = 1.5*max(vecnorm(cgrph.r(:,:)));
x1 = linspace(-L,L,300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];
ntargs = size(targs,2);

% evaluate potentials and return region labels
tic;
[uscat, targdomain] = chunkerkerneval(cgrph, kerns_eval, dens, targs);
uscat = reshape(uscat,size(xx));
tplot = toc

% identify points in computational domain
in = chunkerinterior(cgrph,{x1,x1});
out = ~in;

% get incoming solution
uin = zeros(size(xx));
uin(targdomain==1) = planewave([zks(1);0],targs(:,targdomain==1));

utot = uin + uscat;

umax = max(abs(utot(:))); 
figure(2);clf
h = pcolor(xx,yy,imag(utot)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); clim([-umax,umax]);
hold on 
plot(cgrph,'k')
axis equal
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)





function edge_2_verts = build_loop_edges(id_verts)
edge_2_verts = [id_verts;circshift(id_verts,1)];
end