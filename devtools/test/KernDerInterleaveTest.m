KernDerInterleaveTest0()
KernDerInterleaveTest1()
KernDerInterleaveTest2()

function KernDerInterleaveTest0()
% test the relationship between interleaved Helmholtz kernels

% problem parameters
zk = 1;
coefs = [1;0.5];
coefa = [1,0.3;-1i,0.2];

src = []; src.r = [0;-1]; src.n = [1;-2];
targ = []; targ.r = [[1.1;0.3], [0.3;0.8]]; targ.n = [[-1;0.3],[-1;0.4]];

% kernels
skern = kernel('h','s',zk);
dkern = kernel('h','d',zk);
ckern = kernel('h','c',zk,coefs);

spkern = kernel('h','sp',zk);
dpkern = kernel('h','dp',zk);
cpkern = kernel('h','cp',zk,coefs);

allkern = kernel('h','all',zk,coefa);
trkern = kernel('h','trans_rep',zk,coefs);
trpkern = kernel('h','trans_rep_prime',zk,coefs);
c2trkern = kernel('h','c2tr',zk,coefs);

sgkern = kernel(@(s,t) chnk.helm2d.kern(zk,s,t,'sgrad'));
dgkern = kernel(@(s,t) chnk.helm2d.kern(zk,s,t,'dgrad'));
cgkern = kernel(@(s,t) chnk.helm2d.kern(zk,s,t,'cgrad',coefs));
trgkern = kernel(@(s,t) chnk.helm2d.kern(zk,s,t,'trans_rep_grad',coefs));

% evaluate kernels
sval = skern.eval(src,targ);
dval = dkern.eval(src,targ);
cval = ckern.eval(src,targ);
trval = trkern.eval(src,targ);

spval = spkern.eval(src,targ);
dpval = dpkern.eval(src,targ);
cpval = cpkern.eval(src,targ);
trpval = trpkern.eval(src,targ);

sgval = sgkern.eval(src,targ);
dgval = dgkern.eval(src,targ);
cgval = cgkern.eval(src,targ);
trgval = trgkern.eval(src,targ);

allval = allkern.eval(src,targ);

c2trval = c2trkern.eval(src,targ);

% test combined
assert(norm((coefs(1)*dval + coefs(2)*sval)-cval) <  1e-12)
assert(norm((coefs(1)*dpval + coefs(2)*spval)-cpval) <  1e-12)

% test all
all_assemb = zeros(size(allval));
all_assemb(1:2:end, 1:2:end) = coefa(1,1)*dval;
all_assemb(1:2:end, 2:2:end) = coefa(1,2)*sval;
all_assemb(2:2:end, 1:2:end) = coefa(2,1)*dpval;
all_assemb(2:2:end, 2:2:end) = coefa(2,2)*spval;
assert(norm(all_assemb - allval) <  1e-12)

% test transmission representiatoon
assert(norm([coefs(1)*dval , coefs(2)*sval]-trval) <  1e-12)
assert(norm([coefs(1)*dpval , coefs(2)*spval]-trpval) <  1e-12)

% test gradients
assert(norm((coefs(1)*dgval + coefs(2)*sgval)-cgval) <  1e-12)
assert(norm([coefs(1)*dgval , coefs(2)*sgval]-trgval) <  1e-12)

% compare normal derivate with grad dot n
dpval2 = sum(reshape(targ.n, 2,2).* reshape(dgval,2,[]),1);
assert(norm(dpval - dpval2(:)) <  1e-12)

% test c2trans
c2trval_assemb = zeros(size(c2trval));
c2trval_assemb(1:2:end, :) = coefs(1)*dval + coefs(2)*sval;
c2trval_assemb(2:2:end, :) = coefs(1)*dpval + coefs(2)*spval;
assert(norm(c2trval_assemb-c2trval) <  1e-12)
end

function KernDerInterleaveTest1()
% test the relationship between interleaved Helmholtz difference kernels

% problem parameters
zk = [1,2];
coefs = [[1;0.5],[1;0.5]];
coefa = [1,0.3;-1i,0.2]; coefa = repmat(coefa,[1,1,2]);
coefb = reshape(coefs,2,1,2); coefb = repmat(coefb,[1,2,1]);

src = []; src.r = [0;-1]; src.n = [1;-2];
targ = []; targ.r = [[1.1;0.3], [0.3;0.8]]; targ.n = [[-1;0.3],[-1;0.4]];

% kernels
skern = kernel('helmdiff','s',zk);
dkern = kernel('helmdiff','d',zk);
ckern = kernel('helmdiff','c',zk,coefs);

spkern = kernel('helmdiff','sp',zk);
dpkern = kernel('helmdiff','dp',zk);
cpkern = kernel('helmdiff','cp',zk,coefs);

allkern = kernel('helmdiff','all',zk,coefa);
trkern = kernel('helmdiff','trans_rep',zk,coefs);
trpkern = kernel('helmdiff','trans_rep_prime',zk,coefs);
c2trkern = kernel('helmdiff','c2tr',zk,coefb);

sgkern = kernel(@(s,t) chnk.helm2d.kern(zk(1),s,t,'sgrad_diff') ...
    - chnk.helm2d.kern(zk(2),s,t,'sgrad_diff'));
dgkern = kernel(@(s,t) chnk.helm2d.kern(zk(1),s,t,'dgrad_diff') ...
    - chnk.helm2d.kern(zk(2),s,t,'dgrad_diff'));
cgkern = kernel(@(s,t) chnk.helm2d.kern(zk(1),s,t,'cgrad_diff',coefs) ...
    - chnk.helm2d.kern(zk(2),s,t,'cgrad_diff',coefs));
trgkern = kernel(@(s,t) chnk.helm2d.kern(zk(1),s,t,'trans_rep_grad_diff',coefb(:,:,1)) ...
    - chnk.helm2d.kern(zk(2),s,t,'trans_rep_grad_diff',coefb(:,:,2)));

% evaluate kernels
sval = skern.eval(src,targ);
dval = dkern.eval(src,targ);
cval = ckern.eval(src,targ);
trval = trkern.eval(src,targ);

spval = spkern.eval(src,targ);
dpval = dpkern.eval(src,targ);
cpval = cpkern.eval(src,targ);
trpval = trpkern.eval(src,targ);

sgval = sgkern.eval(src,targ);
dgval = dgkern.eval(src,targ);
cgval = cgkern.eval(src,targ);
trgval = trgkern.eval(src,targ);

allval = allkern.eval(src,targ);

c2trval = c2trkern.eval(src,targ);

% test combined
assert(norm((coefs(1)*dval + coefs(2)*sval)-cval) <  1e-12)
assert(norm((coefs(1)*dpval + coefs(2)*spval)-cpval) <  1e-12)

% test all
all_assemb = zeros(size(allval));
all_assemb(1:2:end, 1:2:end) = coefa(1,1)*dval;
all_assemb(1:2:end, 2:2:end) = coefa(1,2)*sval;
all_assemb(2:2:end, 1:2:end) = coefa(2,1)*dpval;
all_assemb(2:2:end, 2:2:end) = coefa(2,2)*spval;
assert(norm(all_assemb - allval) <  1e-12)

% test transmission representiatoon
assert(norm([coefs(1)*dval , coefs(2)*sval]-trval) <  1e-12)
assert(norm([coefs(1)*dpval , coefs(2)*spval]-trpval) <  1e-12)

% test gradients
assert(norm((coefs(1)*dgval + coefs(2)*sgval)-cgval) <  1e-12)
assert(norm([coefs(1)*dgval , coefs(2)*sgval]-trgval) <  1e-12)

% compare normal derivate with grad dot n
dpval2 = sum(reshape(targ.n, 2,2).* reshape(dgval,2,[]),1);
assert(norm(dpval - dpval2(:)) <  1e-12)

% test c2trans
c2trval_assemb = zeros(size(c2trval));
c2trval_assemb(1:2:end, :) = coefs(1)*dval + coefs(2)*sval;
c2trval_assemb(2:2:end, :) = coefs(1)*dpval + coefs(2)*spval;
assert(norm(c2trval_assemb-c2trval) <  1e-12)
end

function KernDerInterleaveTest2()
% test the relationship between interleaved Laplace kernels

% problem parameters
coefs = [1;0.5];
coefa = [1,0.3;-1i,0.2];

src = []; src.r = [0;-1]; src.n = [1;-2];
targ = []; targ.r = [[1.1;0.3], [0.3;0.8]]; targ.n = [[-1;0.3],[-1;0.4]];

% kernels
skern = kernel('l','s');
dkern = kernel('l','d');
ckern = kernel('l','c',coefs);

spkern = kernel('l','sp');
dpkern = kernel('l','dp');
cpkern = kernel('l','cp',coefs);

sgkern = kernel(@(s,t) chnk.lap2d.kern(s,t,'sgrad'));
dgkern = kernel(@(s,t) chnk.lap2d.kern(s,t,'dgrad'));
cgkern = kernel(@(s,t) chnk.lap2d.kern(s,t,'cgrad',coefs));

% evaluate kernels
sval = skern.eval(src,targ);
dval = dkern.eval(src,targ);
cval = ckern.eval(src,targ);

spval = spkern.eval(src,targ);
dpval = dpkern.eval(src,targ);
cpval = cpkern.eval(src,targ);

sgval = sgkern.eval(src,targ);
dgval = dgkern.eval(src,targ);
cgval = cgkern.eval(src,targ);

% test combined
assert(norm((coefs(1)*dval + coefs(2)*sval)-cval) <  1e-12)
assert(norm((coefs(1)*dpval + coefs(2)*spval)-cpval) <  1e-12)

% test gradients
assert(norm((coefs(1)*dgval + coefs(2)*sgval)-cgval) <  1e-12)

% compare normal derivate with grad dot n
dpval2 = sum(reshape(targ.n, 2,2).* reshape(dgval,2,[]),1);
assert(norm(dpval - dpval2(:)) <  1e-12)

end