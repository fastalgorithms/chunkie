function demo_quasiproxy

% addpath(genpath('/Users/gillmana/Documents/Research/Flatiron/chunkie/'))
% addpath(genpath('/Users/gillmana/Documents/Research/Flatiron/chunkie/+chnk/+flam'))

% Set wave speed in each layer
kh1 = 10;
kh2 = 5;
% angle 
theta = -pi/5;

% test target location for convergence.
xxtrg = [];
xxtrg.r = [0.2; 0.5];
xxtrg.n = [0; 0];

% Number of panels
Npan = 40;

d = 1; % length of unit cell

% top and bottom
ht = 1;
hb = -1;

% look for woods anomaly
K = 20;
% number of terms in all Bragg expansions
KK=2*K+1;
kappa = kh1*cos(theta)+2*pi*[-K:1:K]/d;
ku = sqrt(kh1^2-kappa.^2);
kd = sqrt(kh2^2-kappa.^2);

tol_wood=1e-14;
if min(abs(ku))<tol_wood
    fprintf('\n Warning: Wood anormaly is detected at top layer! (checking tolerance set to %5.2e)',...
        tol_wood);
end
if min(abs(kd))<tol_wood
    fprintf('\n Warning: Wood anormaly is detected at bottom layer! (checking tolerance set to %5.2e)',...
        tol_wood);
end

ima = sqrt(-1);

% Initialize the geometry
cparams = [];
cparams.ta = -0.5;
cparams.tb = 0.5;
chnkr = chunkerfuncuni(@(t) fcurve(t), Npan, cparams);

full_sys = chnk.quasiproxy.build_sys(chnkr,kh1,kh2,d,theta,ht,hb,cparams.ta,K);

% make RHS (only nonzero part)
rhs = chnk.quasiproxy.make_rhs(chnkr,kh1,theta);

% solve the linear system
[proxy_dens,bragg_coef,interface_dens] = chnk.quasiproxy.solve_sys(full_sys,rhs,KK);


%% flux error check
    fprintf('\n All wavespeed are real, energy is preserved: ');
    
    k1x=kh1*cos(theta);
    k1y=kh1*sin(theta);
        [flux_up,flux_down]=chnk.quasiproxy.compute_flux(bragg_coef(1:KK),bragg_coef(KK+1:2*KK),...
        k1x,kh1,kh2,d,K);
    total_flux_up        = sum(real(flux_up));
    total_flux_down      = sum(real(flux_down));
    flux_error           = abs((total_flux_up+total_flux_down-abs(k1y))/abs(k1y));
    fprintf('Flux error est.    = %g\n\n',      flux_error);
    fprintf('\n  theta value , flux error est.    = %g\n',   flux_error );

numlayer = 1;

uapp = chnk.quasiproxy.eval_approx(full_sys,chnkr,interface_dens,proxy_dens,numlayer,d,kh1,theta,xxtrg);    

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
