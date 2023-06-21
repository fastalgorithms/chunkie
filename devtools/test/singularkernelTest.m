%SINGULARKERNELTEST check that principal value and hypersingular type
% quadratures are working.
% 

clearvars; close all;
seed = 8675309;
rng(seed);
addpaths_loc();

doadap = false;

% geometry parameters and construction

cparams = [];
cparams.eps = 1.0e-6;
pref = []; 
pref.k = 10;
narms = 5;
amp = 0.5;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

% sources

ns = 10;
ts = 0.0+2*pi*rand(ns,1);
sources = starfish(ts,narms,amp);
sources = 3.0*sources;
strengths = randn(ns,1);

hold off
plot(chnkr)
hold on
scatter(sources(1,:),sources(2,:),'o')
axis equal 

%

% kernel defs

kernd = kernel('lap','d');
kerns = kernel('lap','s');
kerndgrad = kernel('lap','dgrad');
kernsgrad = kernel('lap','sgrad');

% eval u and dudn on boundary

srcinfo = []; srcinfo.r = sources; 

[ubdry,gradubdry] = kerns.fmm(eps,srcinfo,chnkr.r(:,:),strengths,2);
unbdry = sum(chnkr.n(:,:).*gradubdry,1);

% evaluate gradu on boundary again using singular matrix builders 
% and the identity: 
%      grad u = 2*(grad S dudn - grad D u)
% integral for grad S is interpreted in the principal value sense and the 
% integral for grad D is interpreted in the Hadamard finite part sense

dgradmat = chunkermat(chnkr,kerndgrad);
sgradmat = chunkermat(chnkr,kernsgrad);

gradu = 2*(sgradmat*unbdry(:) - dgradmat*ubdry(:));

relerr = norm(gradu(:)-gradubdry(:))/norm(gradubdry(:));

fprintf('relative frobenius error (Greens ID grad) %5.2e\n',relerr);

assert(relerr < 1e-2);

smat = chunkermat(chnkr,kerns);
mu = smat\ubdry(:);
gradu2 = sgradmat*mu + reshape(chnkr.n(:,:).*mu(:).',2*chnkr.npt,1)*0.5;

relerr = norm(gradu2(:)-gradubdry(:))/norm(gradubdry(:));

fprintf('relative frobenius error (grad S, computed mu) %5.2e\n',relerr);

%%
% need legendre spectral diff matrix to test this part 

[~,~,u,v] = lege.exps(chnkr.k);
dermat  = v(:,1:end-1)*lege.derpol(u);

dmat = chunkermat(chnkr,kernd);
mu = (-0.5*eye(chnkr.npt)+dmat)\ubdry(:);

mu = reshape(mu,chnkr.k,chnkr.nch);
muder = 0.5*(dermat*mu)./chnkr.h(:).';

t = -chnk.perp(chnkr.n(:,:));
gradu3 = dgradmat*mu(:) + reshape(t.*(muder(:).'),2*chnkr.npt,1)*(-0.5);

relerr = norm(gradu3(:)-gradubdry(:))/norm(gradubdry(:));

fprintf('relative frobenius error (grad D, computed mu) %5.2e\n',relerr);

