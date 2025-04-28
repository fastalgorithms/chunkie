complexificationTest0();


function complexificationTest0()
% test solve_sommerfeld_dens and eval_sommerfeld_dens
% and compare to complexification solution

%%%%
%
%       .   .   .   testing a point charge below
%
%%%%

eps = 1E-16;
zks = pi*[1,1.3];
zk1 = zks(1);
zk2 = zks(2);

close all

a = 3;
b = 1/a/2;
t0 =-5;
t1 = 5;

f = @(t) flat_interface(t, a, b, t0, t1);
nch = 50;
xrad = -log(eps)/abs(real(zk1)) + max([abs(t0),abs(t1)]);
xmin = -xrad;
xmax =  xrad;
cparams = [];
cparams.ta = xmin;
cparams.tb = xmax;
cparams.ifclosed = 0;
chnkr = chunkerfuncuni(f, nch, cparams);


ddiff = kernel('helmdiff', 'd', zks);
sdiff = (-1)*kernel('helmdiff', 's', zks);
dpdiff = kernel('helmdiff', 'dp', zks);
spdiff = (-1)*kernel('helmdiff', 'sp', zks);

K = kernel([ddiff, sdiff; ...
           dpdiff, spdiff]);

n = 2*chnkr.npt;
sysmat = chunkermat(chnkr, K);
sysmat = sysmat - eye(n);

%% Analytic solution test

r = chnkr.r;
r0 = [0.1; -3];
s = [];
s.r = r0;

sk1 = kernel('helm', 's', zk1);
dk1 = kernel('helm', 'd', zk1);
skp1 = kernel('helm', 'sp', zk1);
rhs = complex(zeros(n,1));
rhs(1:2:end) = sk1.eval(s, chnkr);
rhs(2:2:end) = skp1.eval(s, chnkr);

dens = sysmat \ rhs;

rt = [0.3; 2];
targ = [];
targ.r = rt;
Keval = kernel([dk1, (-1)*sk1]);

opts = [];
opts.forcesmooth = true;
opts.accel = false;
pot = chunkerkerneval(chnkr, Keval, dens, targ, opts);

uex = sk1.eval(s, targ);
assert(norm(uex - pot) < 1e-10); 





end


function [f,fd,fdd] = flat_interface(t, a, b, t0, t1)
    
    phi   = @(t,u,v,z) u*(t-z).*erfc(u*(t-z))*v - exp(-u^2*(t-z).^2)/sqrt(pi)*v;
    phid  = @(t,u,v,z) u*erfc(u*(t-z))*v;
    phidd = @(t,u,v,z) -u*u*exp(-u^2*(t-z).^2)*2*v/sqrt(pi);
    f = zeros([2,size(t)]);
    fd = zeros([2,size(t)]);
    fdd = zeros([2,size(t)]);
    
    f(1,:) = t + 1i*(phi(t,a,b,t0) - phi(t,-a,b,t1)); 
    fd(1,:)= 1 + 1i*(phid(t,a,b,t0) - phid(t,-a,b,t1));
    fdd(1,:) = 1i*(phidd(t,a,b,t0) - phidd(t,-a,b,t1));
    
    f(2,:) = 0;
    fd(2,:) = 0;
    fdd(2,:) = 0;
        
end

