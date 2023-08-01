%GUIDE_SIMPLEBVPS
%
% This script complements the chunkie guide section on simple boundary 
% problems.
% It shows how to discretize some common PDEs using chunkie
%

addpaths_loc();
rng(1234);

% START SPRIME
kernsp = kernel('lap','sprime');
% END SPRIME

% START LAPLACE NEUMANN
% get a chunker discretization of a starfish domain
chnkr = chunkerfunc(@(t) starfish(t));

% define a boundary condition. because of the symmetries of the 
% starfish, this function integrates to zero
f = @(r) r(1,:).^2.*r(2,:);
rhs = f(chnkr.r); rhs = rhs(:);

% get the kernel 
kernsp = kernel('lap','sprime');

% get a matrix discretization of the boundary integral equation 
sysmat = chunkermat(chnkr,kernsp); % just the sprime part
% add the identity term and "ones matrix"
sysmat = sysmat + 0.5*eye(chnkr.npt) + onesmat(chnkr); 

% solve the system 
sigma = gmres(sysmat,rhs);

% grid for plotting solution
x1 = linspace(-2,2,300);
[xx,yy] = meshgrid(x1,x1);
targs = [xx(:).'; yy(:).'];
in = chunkerinterior(chnkr,targs);
uu = nan(size(xx));

% need the single layer, not it's dervative, to evaluate u
kerns = kernel('lap','s');
uu(in) = chunkerkerneval(chnkr,kerns,sigma,targs(:,in));

% plot
figure(1); clf
h = pcolor(xx,yy,uu); set(h,'EdgeColor','none'); colorbar
% END LAPLACE NEUMANN
saveas(figure(1),"guide_simplebvps_laplaceneumann.png")

