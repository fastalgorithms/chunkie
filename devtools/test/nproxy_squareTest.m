%NPROXY_SQUARE test
%
% Determine number of proxy points using all 3 methods for a helmholtz
% kernel, and a mixed laplace/helmholtz kernel

zk = 3.2;
K = kernel('helm', 's', zk);
width = 1.1;

% Use default method
npxy = chnk.flam.nproxy_square(K, width);

assert (npxy < 50)

% Use manually prescribed wavenumber
opts = [];
opts.npxy_wavenumber = 64;
npxy = chnk.flam.nproxy_square(K, width, opts);

assert (npxy < 200);

% Use id to determine the number of points
opts = [];
opts.npxy_method = 'id';
npxy = chnk.flam.nproxy_square(K, width, opts);

assert(npxy < 50);

% Use integral test to determine the number of points

opts = [];
opts.npxy_method = 'integral';
npxy = chnk.flam.nproxy_square(K, width, opts);

assert(npxy < 150);

% Now try with a composite kernel

zk = 3.2;
K1 = kernel('helm', 's', zk);
K2 = kernel('lap', 'd');
K = kernel([K1, K2]);
width = 1.1;

% Use default method
npxy = chnk.flam.nproxy_square(K, width);

assert (npxy < 50)

% Use manually prescribed wavenumber
opts = [];
opts.npxy_wavenumber = 64;
npxy = chnk.flam.nproxy_square(K, width, opts);

assert (npxy < 200);

% Use id to determine the number of points
opts = [];
opts.npxy_method = 'id';
npxy = chnk.flam.nproxy_square(K, width, opts);

assert(npxy < 50);

% Use integral test to determine the number of points

opts = [];
opts.npxy_method = 'integral';
npxy = chnk.flam.nproxy_square(K, width, opts);

assert(npxy < 150);

%% Now try with a matrix of kernels

zk = 3.2;
zk2 = 1.3;
K1 = kernel('helm', 's', zk);
K2 = kernel('lap', 'd');
K3 = kernel('lap', 'c', [1, 1.3]);
K4 = kernel('helm', 'd', zk2);

K = kernel([K1, K2; K3, K4]);
width = 1.1;

% Use default method
npxy = chnk.flam.nproxy_square(K, width);

assert (npxy < 50)

% Use manually prescribed wavenumber
opts = [];
opts.npxy_wavenumber = 64;
npxy = chnk.flam.nproxy_square(K, width, opts);

assert (npxy < 200);

% Use id to determine the number of points
opts = [];
opts.npxy_method = 'id';
npxy = chnk.flam.nproxy_square(K, width, opts);

assert(npxy < 100);

% Use integral test to determine the number of points

opts = [];
opts.npxy_method = 'integral';
npxy = chnk.flam.nproxy_square(K, width, opts);

assert(npxy < 150);

