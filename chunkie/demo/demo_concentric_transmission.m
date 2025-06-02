%DEMO_CONCENTRIC_TRANSMISSION
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
e2v_out = [id_vert_out;circshift(id_vert_out,1)];

id_vert_in = 5:8;
e2v_in = [id_vert_in;circshift(id_vert_in,1)];

id_vert_in2 = [5,9,7,10];
e2v_in2 =  [id_vert_in2;circshift(id_vert_in2,1)];

id_vert_in3 = [5,11,7,12];
e2v_in3 =  [id_vert_in3;circshift(id_vert_in3,1)];

% combine lists of edges
edge_2_verts = [e2v_out,e2v_in,e2v_in2,e2v_in3];

cparams = [];
cparams.maxchunklen = 4.0/max(zks);
cgrph = chunkgraph(verts,edge_2_verts,[],cparams);
figure(1);
clf
plot_regions(cgrph)
xlim([-2,2])
axis equal 
% END DOMAIN

saveas(figure(1),"demo_concentric_domain.png")

%% Build system
nreg = length(cgrph.regions);

% assign region wavenumbers
ks = [zks(1),repmat(zks(2),1,4),repmat(zks(3),1,2),repmat(zks(4),1,2),zks(5)];

% jump condition coefficients
coefs = ones(1,nreg);


edge_regs = find_edge_regions(cgrph);
% set up parameters for planewave data
opts = [];
opts.bdry_data_type = 'pw';

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
figure(2); clf
h = pcolor(xx,yy,imag(utot)); set(h,'EdgeColor','none'); colorbar
colormap(redblue); clim([-umax,umax]);
hold on 
plot(cgrph,'k')
axis equal
title('$u^{\textrm{tot}}$','Interpreter','latex','FontSize',12)

% END TRANSMISSION PROBLEM
saveas(figure(2),"demo_concentric_fields.png")


