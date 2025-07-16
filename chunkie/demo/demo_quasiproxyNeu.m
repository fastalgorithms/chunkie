function demo_quasiproxyNeu

% addpath(genpath('/Users/gillmana/Documents/Research/Flatiron/chunkie/'))
% addpath(genpath('/Users/gillmana/Documents/Research/Flatiron/chunkie/+chnk/+flam'))

% Set wave speed in each layer
kh1 = 10;
% angle 
theta = -pi/5;

% test target location for convergence.
xxtrg = [];
% xxtrg.r = [0.2; 0.5];
xxtrg.r = [0.2; 1.5];
xxtrg.n = [0; 0];

% Number of panels
Npan = 40;

d = 1; % length of unit cell

% top and bottom
ht = 1;

% look for woods anomaly
K = 20;
% number of terms in all Bragg expansions
KK=2*K+1;
kappa = kh1*cos(theta)+2*pi*[-K:1:K]/d;
ku = sqrt(kh1^2-kappa.^2);

tol_wood=1e-14;
if min(abs(ku))<tol_wood
    fprintf('\n Warning: Wood anormaly is detected at top layer! (checking tolerance set to %5.2e)',...
        tol_wood);
end


% Initialize the geometry
cparams = [];
cparams.ta = -0.5;
cparams.tb = 0.5;
chnkr = chunkerfuncuni(@(t) fcurve(t), Npan, cparams);

full_sys = chnk.quasiproxy.build_sysNeu(chnkr,kh1,d,theta,ht,cparams.ta,K);

% make RHS (only nonzero part)
rhs = chnk.quasiproxy.make_rhsNeu(chnkr,kh1,theta);

% solve the linear system
[proxy_dens,bragg_coef,interface_dens] = chnk.quasiproxy.solve_sysNeu(full_sys,rhs,KK);


uapp = chnk.quasiproxy.eval_approxNeu(full_sys,chnkr,interface_dens,proxy_dens,bragg_coef,KK,d,ht,theta,xxtrg);    

uapp

keyboard


return

% geometry
function r = fcurve(t)
    t1 = t(:);
    r = zeros(2,length(t1));
    r(1,:) = t1.';
    r(2,:) = 0.25*sin(2*pi*(t1+0.5));
    r = reshape(r, [2, size(t)]);
return
