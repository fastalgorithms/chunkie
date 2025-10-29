gmresF2KconvergenceTest0();


function gmresF2KconvergenceTest0()

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple chunker 
%            and tests new gmres using chunkermat
%            from the Helmholtz transmission problem
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iseed = 8675309;
rng(iseed);

kvec = 20*[1;-1.5];
zk = norm(kvec);

% discretize domain
cparams = [];
cparams.eps = 1.0e-6;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zk; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.25;
chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
chnkr = refine(chnkr);

fkern = kernel('helm','c',zk,[1,-zk*1i]);
sysmatK = chunkermat(chnkr,fkern)*1e-6;
sysmatA = 0.5*eye(chnkr.k*chnkr.nch) + sysmatK;
rhs = -planewave(kvec(:),chnkr.r(:,:)); rhs = rhs(:);

tol = 1e-40;
maxit = 200;
opts.tol = tol;
opts.maxit = maxit;
dens0 = sysmatA\rhs;
[dens1,~,~,iter1,resvec1] = gmres(sysmatA,rhs,[],tol,maxit); 
[dens2,iter2,resvec2,~] = gmresF2K(sysmatK,rhs,0.5,[],opts); 

clf;
plot(log10(abs(resvec1)));
hold on; 
plot(log10(abs(resvec2)));
legend('GMRES','GMRESF2K');
title('convergence of gmres vs gmresF2K');
norm(dens0-dens1)
norm(dens0-dens2)
assert(norm(dens0-dens1)<1e-10);
assert(norm(dens0-dens2)<1e-10);

residue1 = sysmatA*dens1-rhs;
residue2 = sysmatA*dens2-rhs;

assert(norm(residue1)<1e-10);
assert(norm(residue2)<1e-10);

fprintf('standard gmres converged in %d iterations to %.4e.\n',iter1(end),resvec1(end));
fprintf('gmres for zI+K converged in %d iterations to %.4e.\n',iter2,resvec2(end));

end
