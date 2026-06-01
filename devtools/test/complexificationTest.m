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

a = 3;
b = 1;
t0 =-5;
t1 = 5;

f = @(t) chnk.curves.complex_interface(t, a, b, t0, t1);
nch = 30;
xrad = -log(eps)/abs(real(zk1)) + max([abs(t0),abs(t1)]);
xmin = -xrad;
xmax =  xrad;
cparams = [];
cparams.ta = xmin;
cparams.tb = xmax;
cparams.ifclosed = 0;
chnkr = chunkerfuncuni(f, nch, cparams);

rends = chunkends(chnkr,[1,chnkr.nch]);
verts = rends(:,[1,4]);
fchnk{1} = chnkr;
cgrph = chunkgraph(real(verts),[1;2],fchnk);

Kres = kernel('lap', 's');
srcinfo = [];
srcinfo.r = [0;0.1];
opts_refine = [];
opts_refine.targfun = @(t) Kres.eval(srcinfo, t);
lastwarn('','');
cgrph = refine(cgrph, opts_refine);
[~, warnID] = lastwarn();
assert(isempty(warnID), 'Warning in complex chunkgraph with refine');


ddiff = kernel('helmdiff', 'd', zks);
sdiff = (-1)*kernel('helmdiff', 's', zks);
dpdiff = kernel('helmdiff', 'dp', zks);
spdiff = (-1)*kernel('helmdiff', 'sp', zks);

K = kernel([ddiff, sdiff; ...
           dpdiff, spdiff]);

n = 2*cgrph.npt;
sysmat = chunkermat(cgrph, K);
sysmat = sysmat - eye(n);

%% Analytic solution test

r = cgrph.r;
r0 = [0.1; -3];
s = [];
s.r = r0;

sk1 = kernel('helm', 's', zk1);
dk1 = kernel('helm', 'd', zk1);
skp1 = kernel('helm', 'sp', zk1);
rhs = complex(zeros(n,1));
rhs(1:2:end) = sk1.eval(s, cgrph);
rhs(2:2:end) = skp1.eval(s, cgrph);

dens = sysmat \ rhs;

rt = [[0.3; 2] [-0.5;0.001]];
targ = [];
targ.r = rt;
Keval = kernel([dk1, (-1)*sk1]);

opts = [];
% opts.forcesmooth = true;
opts.accel = false;
pot = chunkerkerneval(cgrph, Keval, dens, targ, opts);

uex = sk1.eval(s, targ);
assert(norm(uex - pot) < 1e-9); 

end
