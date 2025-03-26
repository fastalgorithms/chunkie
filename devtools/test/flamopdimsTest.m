
zk1 = 10;   % exterior wave number
zk2 = 15;   % coating wave number


zkuse = max(real([zk1, zk2]));
cparams = [];
cparams.eps = 1.0e-10;
cparams.nover = 0;
cparams.maxchunklen = 4.0/zkuse; % setting a chunk length helps when the
                              % frequency is known
                              
pref = []; 
pref.k = 16;
narms = 5;
amp = 0.15;

chnkr1 = chunkerfunc(@(t) starfish(t,narms,amp, [], [], 1.15),cparams,pref); 

narms2 = 10;
amp2 = 0.03;
chnkr2 = chunkerfunc(@(t) starfish(t, narms2, amp2, [], [], 0.9), cparams, pref); 


ifflam = true;
opts = [];
opts.ifflam = ifflam;

%% Test analytic solution
tic, A = get_linop(chnkr1, chnkr2, [zk1, zk2], opts); toc

src1 = zeros(2,2);
src1(1,1) = 0.01;
src1(2,1) = 0.03;

src1(1,2) = -1;
src1(2,2) = -0.5;

src2 = zeros(2,2);
src2(:,1) = src1(:,1);
src2(1,2) = 5.3;
src2(2,2) = -3.1;

err1  = test_analytic_soln([zk1, zk2], chnkr1, chnkr2, src1, src2, A, ifflam);
assert(err1 < 1e-8);



function err1 = test_analytic_soln(zks, chnkr1, chnkr2, src1 ,src2, A, ifflam)
    zk1 = zks(1);
    zk2 = zks(2);

    dk1 = kernel('helm', 'd', zk1);
    sk1 = kernel('helm', 's', zk1);
    sk1p = kernel('helm', 'sp', zk1);

    dk2 = kernel('helm', 'd', zk2);
    sk2 = kernel('helm', 's', zk2);
    sk2p = kernel('helm', 'sp', zk2);
    ck2 = kernel('helm', 'c', zk2, [2, -2*1j*zk2]);

    [~, n1] = size(src1);
    strengths1 = randn(n1,1) + 1j*randn(n1,1);

    srcinfo1 = [];
    srcinfo1.r = src1;

    targinfo1 = [];
    targinfo1.r = chnkr1.r(:,:);
    targinfo1.n = chnkr1.n(:,:);

    targinfo2 = [];
    targinfo2.r = chnkr2.r(:,:);
    targinfo2.n = chnkr2.n(:,:);

    u1 = sk1.eval(srcinfo1, targinfo1)*strengths1;
    dudn1 = sk1p.eval(srcinfo1, targinfo1)*strengths1;


    srcinfo2 = [];
    srcinfo2.r = src2;
    [~, n2] = size(src2);

    strengths2 = randn(n2,1) + 1j*randn(n2,1);

    u2_1 = sk2.eval(srcinfo2, targinfo1)*strengths2;
    dudn2_1 = sk2p.eval(srcinfo2, targinfo1)*strengths2;

    u2_2 = sk2.eval(srcinfo2, targinfo2)*strengths2;

    ntot = chnkr1.npt*2 + chnkr2.npt;
    n = 2*chnkr1.npt;
    rhs = zeros(ntot,1);
    rhs(1:2:n) = u1 - u2_1;
    rhs(2:2:n) = dudn1 - dudn2_1;
    rhs((n+1):end) = u2_2;

    if ~ifflam
        sol = A\rhs;
    else
        sol = rskelf_sv(A, rhs);
    end

% Now test solution at interior and exterior point

    targ1 = src2(:,2);
    targ2 = src1(:,2);
    targinfo1 = [];
    targinfo1.r = targ1;

    targinfo2 = [];
    targinfo2.r = targ2;

    uex1 = sk1.eval(srcinfo1, targinfo1)*strengths1;
    uex2 = sk2.eval(srcinfo2, targinfo2)*strengths2;


    K1 = kernel([dk1, -1*sk1]);
    u1 = chunkerkerneval(chnkr1, K1, sol(1:n), targinfo1);

    K2 = kernel([dk2, -1*sk2]);
    u2 = chunkerkerneval(chnkr1, K2, sol(1:n), targinfo2);
    u2 = u2 + chunkerkerneval(chnkr2, ck2, sol((n+1):end), targinfo2);
    err1 = norm(u1 - uex1) + norm(u2 - uex2);

    fprintf('error in exterior=%d\n', norm(u1-uex1));
    fprintf('error in coating=%d\n', norm(u2-uex2));
end

function [A] = get_linop(chnkr1, chnkr2, zks, opts)
    zk1 = zks(1);
    zk2 = zks(2);
    if nargin < 4
        opts = [];
    end
    ifflam = false;
    if isfield(opts, 'ifflam')
        ifflam = opts.ifflam;
    end
    chnkrs(2,1) = chunker();
    chnkrs(1,1) = chnkr1;
    chnkrs(2,1) = chnkr2;

%% Define kernels
    skdiff = kernel('helmdiff', 's', [zk1, zk2]);
    skpdiff = kernel('helmdiff', 'sp', [zk1, zk2]); 
    dkdiff = kernel('helmdiff', 'd', [zk1, zk2]);
    dkpdiff = kernel('helmdiff', 'dp', [zk1, zk2]);

    K = kernel([dkdiff, -1*skdiff; ...
            dkpdiff, -1*skpdiff]);


    dk2 = kernel('helm', 'd', zk2);
    sk2 = kernel('helm', 's', zk2);
 
    ck2 = kernel('helm', 'c', zk2, [2, -2*1j*zk2]);
    ck2p = kernel('helm', 'cp', zk2, [2, -2*1j*zk2]);

    K2 = kernel([dk2, -1*sk2]);
    K3 = -1*kernel([ck2; ck2p]);

    Kmat(2,2) = kernel();
    Kmat(1,1) = K;
    Kmat(2,1) = K2;
    Kmat(2,2) = ck2;
    Kmat(1,2) = K3;

%% Evaluate matrix

opts_loc = [];
opts_loc.adaptive_correction = true;

if ~ifflam
    A = chunkermat(chnkrs, Kmat, opts_loc);
    ntot = size(A);
    A = A + eye(ntot);
else
    A = chunkerflam(chnkrs, Kmat, 1, opts_loc); 
end

end
