
zk0 = 2.0; % exterior
zk1 = 3.0; % interior
zk  = [zk1,zk0];

kernparam.zk0 = zk0; kernparam.zk1 = zk1;

cparams = [];   
cparams.nover = 0;  
cparams.ifclosed = 1; 
cparams.eps = 1.0e-10;
modes = 1;

nch = 20;
chnkr = chunkerfuncuni(@(t) chnk.curves.bymode(t,modes,[0,0],[1.2,1]),nch,cparams); 
chnkr = sort(chnkr);



ns = 1;
ts0 = 2*pi*rand(ns,1);
ts1 = 2*pi*rand(ns,1);
src_ext = [0;0];
src_int = [-10;5];

rmin = min(chnkr); rmax = max(chnkr);
xl = rmax(1)-rmin(1); yl = rmax(2)-rmin(2);
nplot = 100;
xtarg = linspace(rmin(1)-xl,rmax(1)+xl,nplot); 
ytarg = linspace(rmin(2)-yl,rmax(2)+yl,nplot);
[xxtarg,yytarg] = meshgrid(xtarg,ytarg);
targets = zeros(2,length(xxtarg(:)));
targets(1,:) = xxtarg(:); targets(2,:) = yytarg(:);
targetsinfo = []; targetsinfo.r = targets;

in = targets(1,:).^2 + targets(2,:).^2 < 0.81; 
out = targets(1,:).^2 + targets(2,:).^2 > 1.21; 


srcinfo0 = []; srcinfo0.r = src_ext;
srcinfo1 = []; srcinfo1.r = src_int;


f_u0 = @(t) chnk.helm2d.kern(zk0,srcinfo0,t,'s');
f_u0n = @(t) chnk.helm2d.kern(zk0,srcinfo0,t,'sprime');
f_u1 = @(t) chnk.helm2d.kern(zk1,srcinfo1,t,'s');
f_u1n = @(t) chnk.helm2d.kern(zk1,srcinfo1,t,'sprime');

u0bdr = f_u0(chnkr);
u0nbdr = f_u0n(chnkr);

u1bdr = f_u1(chnkr);
u1nbdr = f_u1n(chnkr);

% the incident field at the boundary,

wts = chnkr.wts;
wts = wts(:);
swts = sqrt(wts);
npt = chnkr.npt;
nn = 2*npt;
rhs = complex(zeros(2*npt,1));
u_inc =  -u0bdr + u1bdr;
un_inc = -u0nbdr + u1nbdr;


rhs(1:2:nn) = -u_inc.*swts;
rhs(2:2:nn) = -un_inc.*swts/(1j*zk0+1j*zk1)*2;


% % constructing the matrix. 

% D0 - D1    |     ik1S1 - ik0S0
% D0' - D1'  |     ik1S1' - ik0S0' /(1j*zk0+1j*zk1)*2;


fkern_d_diff = @(s,t) chnk.helm2d.kern(zk0,s,t,'D')...
                    - chnk.helm2d.kern(zk1,s,t,'D');
fkern_s_diff = @(s,t) -1j*zk0*chnk.helm2d.kern(zk0,s,t,'S')...
                    + 1j*zk1*chnk.helm2d.kern(zk1,s,t,'S');
fkern_dprime_diff = @(s,t) chnk.helm2d.kern(zk0,s,t,'dprime')...
                         - chnk.helm2d.kern(zk1,s,t,'dprime');
kern_sprime_diff = @(s,t) -1j*zk0*chnk.helm2d.kern(zk0,s,t,'sprime') ...
                        + 1j*zk1*chnk.helm2d.kern(zk1,s,t,'sprime');
% % Construct manually l2scaled matrices

A = zeros(2*npt,'like', 1j);


A11 = chunkermat(chnkr,fkern_d_diff) +  eye(npt);
A12 = chunkermat(chnkr,fkern_s_diff);
A21 = chunkermat(chnkr,fkern_dprime_diff)/(1j*zk0+1j*zk1)*2;
A22 = chunkermat(chnkr,kern_sprime_diff)/(1j*zk0+1j*zk1)*2 + eye(npt);

dd = diag(swts);
ddinv = diag(1.0./swts);
A(1:2:nn,1:2:nn) = dd*A11*ddinv;
A(1:2:nn,2:2:nn) = dd*A12*ddinv;
A(2:2:nn,1:2:nn) = dd*A21*ddinv;
A(2:2:nn,2:2:nn) = dd*A22*ddinv;

% % Use l2scale option of chunkermat
opts = [];
opts.l2scale = 'true';
Ac = zeros(2*npt,'like', 1j);
A11 = chunkermat(chnkr,fkern_d_diff,opts) +  eye(npt);
A12 = chunkermat(chnkr,fkern_s_diff,opts);
A21 = chunkermat(chnkr,fkern_dprime_diff,opts)/(1j*zk0+1j*zk1)*2;
A22 = chunkermat(chnkr,kern_sprime_diff,opts)/(1j*zk0+1j*zk1)*2 + eye(npt);
Ac(1:2:nn,1:2:nn) = A11;
Ac(1:2:nn,2:2:nn) = A12;
Ac(2:2:nn,1:2:nn) = A21;
Ac(2:2:nn,2:2:nn) = A22;



% % Compare the two matrices
x = A\rhs;
x2 = Ac\rhs;

err1 = norm(x-x2);
err2 = norm(A-Ac);
fprintf('error in matrix norm=%d\n',err2)
fprintf('error in density=%d\n',err1);


