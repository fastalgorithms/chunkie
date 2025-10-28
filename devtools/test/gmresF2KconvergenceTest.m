clearvars; close all;
iseed = 8675309;
rng(iseed);
npts = 1e3;
% A = zI+K
matK = (rand(npts)-0.5)*1e-8; % K is a small random matrix
z = 1i; % z is any complex number
matA = z*eye(npts) + matK;
rhs = matA * rand(npts,1);
rhs = rhs(:);

start = tic; [sol1,flag,relres,iter1,resvec1] = gmres(matA,rhs,[],1e-15,100); t1 = toc(start);
% test new gmres 
% set up parameters
opts = [];
opts.max_iterations = 1e4;
opts.threshold = 1e-100;
x = zeros(size(rhs));
% test gmres for 
[sol3,iter3,resvec3] = gmresF2K(z,matK,rhs,x, opts);
resvec3 = abs(resvec3);
resvec3 = resvec3(resvec3~=0);

clf;
plot(log10(resvec1),Color='b',LineWidth=2); 
hold on;
plot(log10(resvec3),Color='r',LineWidth=2); 
legend('Native GMRES', 'GMRES for F2K operators');

fprintf("Native GMRES exited with flag %i (3=stalled), converged to %f digits. \n",flag,log10(resvec1(end)));
fprintf("New GMRES for F2K operators converged to %f digits. \n",log10(resvec3(end)));
assert(resvec3(end)<opts.threshold)
