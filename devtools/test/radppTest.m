%RADPPTEST test kernel.radgreen and kernel.radpp

rng(42);

test_fit_limits();
test_tensors();
test_formats();
test_laplace();
test_helmholtz();
test_flex();
test_dim3();
test_radpp();
test_sing();
test_arithmetic();
test_options();
test_edge();
test_radpp_formats();
test_warnings();
test_chunkermat();
test_bad_inputs();


function test_fit_limits()
% nlevmax warns and keeps the panels, maxsub errors

f = @(r) 1./r;
lastwarn('');
warning('off','CHNK:radpp:fit:nlevmax');
info = chnk.radpp.fit(f,1e-6,1,struct('nlevmax',3));
warning('on','CHNK:radpp:fit:nlevmax');
[~,id] = lastwarn();
assert(strcmp(id,'CHNK:radpp:fit:nlevmax'),'radpp fit: no nlevmax warning');
assert(info.nbin == 4,'radpp fit: unresolved panels not kept');

caught = false;
try
    chnk.radpp.fit(f,1e-6,1,struct('maxsub',10));
catch
    caught = true;
end
assert(caught,'radpp fit: maxsub not an error');

end


function test_tensors()
% tensors up to fourth order from the on-axis profiles of a gaussian

for dim = [2 3]
    green = @(s,t) gaussgreen(dim,s,t);
    n = 300;
    X = 1.5*randn(n,3);
    r = sqrt(sum(X(:,1:dim).^2,2));

    out = cell(1,5);
    [out{:}] = green(zeros(dim,1),X(:,1:dim).');
    C = out{1}(:);
    for k = 1:4
        C = [C, chnk.radpp.unpack(out{k+1},dim,k)];
    end

    prof = chnk.radpp.profiles(green,r,4,dim);
    P = chnk.radpp.rotate(X(:,1),X(:,2),(dim==3)*X(:,3),prof,4,dim);
    assert(isequal(size(P),size(C)),'radpp tensors %dd: wrong size',dim);
    err = max(abs(P(:) - C(:)));
    assert(err < 1e-12,'radpp tensors %dd: %.2e',dim,err);
end

end


function test_formats()
% packed and full derivative arrays give the same kernels

[src,targ] = geom();
nu = 0.3;
full = @(s,t) gaussgreen(2,s,t);
types = {'s','sp','dp','clamped_bc','clamped_eval', ...
    'supported_bc','free_bc'};
for i = 1:numel(types)
    a = {};
    if any(strcmp(types{i},{'supported_bc','free_bc'})), a = {nu}; end
    A = kernel.radgreen(full,types{i},a{:}).eval(src,targ);
    B = kernel.radgreen(@(s,t) gausspacked(2,s,t),types{i},a{:}).eval(src,targ);
    err = norm(A - B,'fro')/norm(A,'fro');
    assert(err < 1e-14,'radgreen packed format %s: %.2e',types{i},err);
end

end


function test_laplace()
% layer potentials against chnk.lap2d.kern

[src,targ] = geom();
green = @(s,t) chnk.lap2d.green(s,t);
cf = [0.7 1.3];

tab = {'s','s',0; 'd','d',0; 'sp','sprime',0; 'dp','dprime',0; ...
    'c','c',1; 'cp','cprime',1; 'sg','sgrad',0; 'dg','dgrad',0; ...
    'cg','cgrad',1; 'st','stau',0; 'hilb','hilb',0};
for i = 1:size(tab,1)
    a = {};
    if tab{i,3}, a = {cf}; end
    K = kernel.radgreen(green,tab{i,1},a{:});
    A = K.eval(src,targ);
    B = chnk.lap2d.kern(src,targ,tab{i,2},a{:});
    err = norm(A - B,'fro')/norm(B,'fro');
    assert(isequal(size(A),size(B)),'radgreen lap %s: wrong size',tab{i,1});
    assert(err < 1e-11,'radgreen lap %s: %.2e',tab{i,1},err);
end

% the two row kernels are the rows of s, sp and d, dp
S = kernel.radgreen(green,'s').eval(src,targ);
Sp = kernel.radgreen(green,'sp').eval(src,targ);
D = kernel.radgreen(green,'d').eval(src,targ);
Dp = kernel.radgreen(green,'dp').eval(src,targ);
A = kernel.radgreen(green,'s2trans').eval(src,targ);
assert(norm(A(1:2:end,:) - S,'fro') + norm(A(2:2:end,:) - Sp,'fro') < 1e-12, ...
    'radgreen s2trans');
A = kernel.radgreen(green,'d2trans').eval(src,targ);
assert(norm(A(1:2:end,:) - D,'fro') + norm(A(2:2:end,:) - Dp,'fro') < 1e-12, ...
    'radgreen d2trans');


% custom kernels
cfun = @(D,s,t) cf(1)*D({t.n},{s.n}) + cf(2)*D({t.n},{});
A = kernel.radgreen(green,cfun,'nderiv',2).eval(src,targ);
B = kernel.radgreen(green,'cp',cf).eval(src,targ);
assert(norm(A - B,'fro') < 1e-12*norm(B,'fro'),'radgreen custom scalar');
sfun = @(D,s,t) cat(3,D({},{}),D({t.n},{}));
K = kernel.radgreen(green,sfun,'nderiv',1,'opdims',[2 1]);
A = K.eval(src,targ);
B = kernel.radgreen(green,'s2trans').eval(src,targ);
assert(norm(A - B,'fro') < 1e-12*norm(B,'fro'),'radgreen custom vector');
assert(strcmp(K.sing,'log'),'radgreen custom sing');

end


function test_helmholtz()
% layer potentials and transmission systems against chnk.helm2d.kern

[src,targ] = geom();
zk = 1.1;
green = @(s,t) chnk.helm2d.green(zk,s,t);
cf = [0.7 1.3i];
cf22 = [1 2i; 3 4];

tab = {'s','s',0; 'd','d',0; 'sp','sprime',0; 'dp','dprime',0; ...
    'c','c',1; 'cp','cprime',1; 'sg','sgrad',0; 'dg','dgrad',0; ...
    'cg','cgrad',1; 'c2trans','c2trans',2; 'trans','all',2; ...
    'trans_rep','trans_rep',1; 'trans_rep_prime','trans_rep_prime',1; ...
    'trans_rep_grad','trans_rep_grad',1};
for i = 1:size(tab,1)
    a = {};
    if tab{i,3} == 1, a = {cf}; end
    if tab{i,3} == 2, a = {cf22}; end
    K = kernel.radgreen(green,tab{i,1},a{:});
    A = K.eval(src,targ);
    B = chnk.helm2d.kern(zk,src,targ,tab{i,2},a{:});
    err = norm(A - B,'fro')/norm(B,'fro');
    assert(isequal(size(A),size(B)),'radgreen helm %s: wrong size',tab{i,1});
    assert(err < 1e-11,'radgreen helm %s: %.2e',tab{i,1},err);
    assert(isequal(K.opdims,[size(A,1)/numel(targ.n(1,:)), ...
        size(A,2)/numel(src.n(1,:))]),'radgreen helm %s: opdims',tab{i,1});
end

end


function test_flex()
% plate kernels against chnk.flex2d.kern

[src,targ] = geom();
zk = 1.3;
nu = 0.3;
green = @(s,t) flexgreen(zk,s,t);

tab = {'clamped_bc','clamped_plate_bcs',0; ...
    'clamped_eval','clamped_plate_eval',0; ...
    'supported_bc','supported_plate_bcs',1; ...
    'supported_eval','supported_plate_eval',1; ...
    'free_bc','free_plate_bcs',1; 'free_eval','free_plate_eval',1};
for i = 1:size(tab,1)
    a = {};
    if tab{i,3}, a = {nu}; end
    A = kernel.radgreen(green,tab{i,1},a{:}).eval(src,targ);
    B = chnk.flex2d.kern(zk,src,targ,tab{i,2},a{:});
    err = norm(A - B,'fro')/norm(B,'fro');
    assert(isequal(size(A),size(B)),'radgreen flex %s: wrong size',tab{i,1});
    assert(err < 1e-9,'radgreen flex %s: %.2e',tab{i,1},err);
end

end


function test_dim3()
% three dimensional kernels against contractions of the full tensors

ns = 7;  nt = 9;
src.r = randn(3,ns);
src.n = randn(3,ns);  src.n = src.n./vecnorm(src.n);
targ.r = 3 + randn(3,nt);
targ.n = randn(3,nt);  targ.n = targ.n./vecnorm(targ.n);
green = @(s,t) gaussgreen(3,s,t);
cf = [0.7 1.3i];
cf22 = [1 2i; 3 4];

[v,g,h] = green(src.r,targ.r);
SN = repmat(permute(src.n,[3 2 1]),[nt 1 1]);
TN = repmat(permute(targ.n,[2 3 1]),[1 ns 1]);
D = -tcon(g,{SN});
Sp = tcon(g,{TN});
Dp = -tcon(h,{TN,SN});
sg = zeros(3*nt,ns);
dg = zeros(3*nt,ns);
for i = 1:3
    sg(i:3:end,:) = g(:,:,i);
    for j = 1:3
        dg(i:3:end,:) = dg(i:3:end,:) - h(:,:,i,j).*SN(:,:,j);
    end
end

ref = struct('s',v,'d',D,'sp',Sp,'dp',Dp, ...
    'c',cf(1)*D + cf(2)*v,'cp',cf(1)*Dp + cf(2)*Sp, ...
    'sg',sg,'dg',dg,'cg',cf(1)*dg + cf(2)*sg);

ref.s2trans = zeros(2*nt,ns);
ref.s2trans(1:2:end,:) = v;
ref.s2trans(2:2:end,:) = Sp;
ref.d2trans = zeros(2*nt,ns);
ref.d2trans(1:2:end,:) = D;
ref.d2trans(2:2:end,:) = Dp;
ref.c2trans = zeros(2*nt,ns);
ref.c2trans(1:2:end,:) = cf22(1,1)*D + cf22(1,2)*v;
ref.c2trans(2:2:end,:) = cf22(2,1)*Dp + cf22(2,2)*Sp;
ref.trans = zeros(2*nt,2*ns);
ref.trans(1:2:end,1:2:end) = cf22(1,1)*D;
ref.trans(1:2:end,2:2:end) = cf22(1,2)*v;
ref.trans(2:2:end,1:2:end) = cf22(2,1)*Dp;
ref.trans(2:2:end,2:2:end) = cf22(2,2)*Sp;
ref.trans_rep = zeros(nt,2*ns);
ref.trans_rep(:,1:2:end) = cf(1)*D;
ref.trans_rep(:,2:2:end) = cf(2)*v;
ref.trans_rep_prime = zeros(nt,2*ns);
ref.trans_rep_prime(:,1:2:end) = cf(1)*Dp;
ref.trans_rep_prime(:,2:2:end) = cf(2)*Sp;
ref.trans_rep_grad = zeros(3*nt,2*ns);
ref.trans_rep_grad(:,1:2:end) = cf(1)*dg;
ref.trans_rep_grad(:,2:2:end) = cf(2)*sg;

types = fieldnames(ref);
for i = 1:numel(types)
    a = targs(types{i});
    K = kernel.radgreen(green,types{i},a{:},'dim',3);
    A = K.eval(src,targ);
    err = norm(A - ref.(types{i}),'fro')/norm(ref.(types{i}),'fro');
    assert(isequal(size(A),size(ref.(types{i}))),'radgreen 3d %s: size',types{i});
    assert(err < 1e-12,'radgreen 3d %s: %.2e',types{i},err);
end

% packed derivative arrays
packed = @(s,t) gausspacked(3,s,t);
for type = {'dp','cg','trans'}
    a = targs(type{1});
    A = kernel.radgreen(green,type{1},a{:},'dim',3).eval(src,targ);
    B = kernel.radgreen(packed,type{1},a{:},'dim',3).eval(src,targ);
    assert(norm(A - B,'fro') < 1e-14*norm(A,'fro'),'radgreen 3d packed %s',type{1});
end

end


function test_radpp()
% the piecewise polynomial fit against exact evaluation

[src,targ] = geom();
nu = 0.3;
zk = 1.1;
types = {'s','sp','dp','clamped_bc','clamped_eval', ...
    'supported_bc','supported_eval','free_bc','free_eval'};
greens = {@(s,t) gaussgreen(2,s,t), @(s,t) chnk.helm2d.green(zk,s,t)};
for j = 1:numel(greens)
    for i = 1:numel(types)
        if j == 2 && any(strcmp(types{i},{'clamped_eval', ...
                'supported_eval','free_bc','free_eval'}))
            continue
        end
        a = {};
        if any(strcmp(types{i},{'supported_bc','supported_eval', ...
                'free_bc','free_eval'})), a = {nu}; end
        Kr = kernel.radpp(greens{j},[1e-6 6],types{i},a{:});
        Kg = kernel.radgreen(greens{j},types{i},a{:});
        A = Kr.eval(src,targ);
        B = Kg.eval(src,targ);
        err = norm(A - B,'fro')/norm(B,'fro');
        assert(err < 1e-10,'radpp %s (green %d): %.2e',types{i},j,err);
        assert(strcmp(Kr.sing,Kg.sing),'radpp %s: sing',types{i});
    end
end

% three dimensions
ns = 7;  nt = 9;
s3.r = 0.5*randn(3,ns);
s3.n = randn(3,ns);  s3.n = s3.n./vecnorm(s3.n);
t3.r = 0.5 + 0.5*randn(3,nt);
t3.n = randn(3,nt);  t3.n = t3.n./vecnorm(t3.n);
green = @(s,t) gaussgreen(3,s,t);
for type = {'s','d','sp','dp','sg','dg','cg','c2trans','trans'}
    a = targs(type{1});
    Kr = kernel.radpp(green,[1e-6 6],type{1},a{:},'dim',3);
    Kg = kernel.radgreen(green,type{1},a{:},'dim',3);
    B = Kg.eval(s3,t3);
    err = norm(Kr.eval(s3,t3) - B,'fro')/norm(B,'fro');
    assert(err < 1e-10,'radpp 3d %s: %.2e',type{1},err);
end

% NaN outside the fit
K = kernel.radpp(greens{1},[1e-3 6],'s');
A = K.eval(struct('r',[10;0]),struct('r',[0;0]));
assert(isnan(A),'radpp: no NaN outside the fit');

% the sum of two greens
g1 = @(s,t) chnk.helm2d.green(zk,s,t);
g2 = @(s,t) chnk.lap2d.green(s,t);
gs = chnk.radpp.sumgreens([2 -1],g1,g2);
A = kernel.radgreen(gs,'sp').eval(src,targ);
B = 2*kernel.radgreen(g1,'sp').eval(src,targ) - ...
    kernel.radgreen(g2,'sp').eval(src,targ);
assert(norm(A - B,'fro') < 1e-12*norm(B,'fro'),'radgreen: sumgreens');

end


function test_sing()
% singularity types of the kernels

green = @(s,t) chnk.lap2d.green(s,t);
tab = {'s','log'; 'sp','log'; 'd','log'; 'dp','hs'; 'sg','pv'; 'dg','hs'; ...
    'st','pv'; 'hilb','pv'};
for i = 1:size(tab,1)
    K = kernel.radgreen(green,tab{i,1});
    assert(strcmp(K.sing,tab{i,2}),'radgreen sing %s: %s',tab{i,1},K.sing);
end

green = @(s,t) chnk.helm2d.green(1.1,s,t);
tab = {'s','log'; 'sp','log'; 'd','log'; 'dp','hs'; 'c','log'; 'cp','hs'; ...
    'sg','pv'; 'dg','hs'; 'cg','hs'; 'c2trans','hs'; 'trans','hs'; ...
    'trans_rep','log'; 'trans_rep_prime','hs'; 'trans_rep_grad','hs'};
for i = 1:size(tab,1)
    a = targs(tab{i,1});
    K = kernel.radgreen(green,tab{i,1},a{:});
    assert(strcmp(K.sing,tab{i,2}),'radgreen helm sing %s: %s',tab{i,1},K.sing);
end

K = kernel.radgreen(green,'s','sing','hs');
assert(strcmp(K.sing,'hs'),'radgreen: sing option');

green = @(s,t) flexgreen(1.3,s,t);
K = kernel.radgreen(green,'clamped_bc');
assert(strcmp(K.sing,'log'),'radgreen sing clamped_bc: %s',K.sing);

end


function test_arithmetic()
% kernel arithmetic, interleaving and dispatch with radgreen and radpp

[src,targ] = geom();
green = @(s,t) chnk.lap2d.green(s,t);
Ks = kernel.radgreen(green,'s');
Kd = kernel.radgreen(green,'d');
Kdp = kernel.radgreen(green,'dp');
Kr = kernel.radpp(green,[1e-6 6],'sp');
Kb = kernel('lap','d');
S = Ks.eval(src,targ);
D = Kd.eval(src,targ);
Sp = Kr.eval(src,targ);
a = 2 + 1i;

chk = @(K,B) norm(K.eval(src,targ) - B,'fro') < 1e-10*norm(B,'fro');
assert(chk(Ks + Kd,S + D),'radgreen plus');
assert(chk(Ks - Kd,S - D),'radgreen minus');
assert(chk(-Ks,-S),'radgreen uminus');
assert(chk(a.*Ks,a*S),'radgreen times');
assert(chk(Ks/a,S/a),'radgreen divide');
assert(chk(Ks + Kb,S + Kb.eval(src,targ)),'radgreen plus built-in');
assert(chk(Ks + Kr,S + Sp),'radgreen plus radpp');
Ksum = Ks + Kd;
assert(strcmp(Ksum.sing,'log'),'radgreen plus sing');
Ksum = Ks + Kdp;
assert(strcmp(Ksum.sing,'hs'),'radgreen plus sing hs');

K = kernel([Ks; Kd]);
A = K.eval(src,targ);
assert(isequal(K.opdims,[2 1]),'radgreen interleave opdims');
assert(norm(A(1:2:end,:) - S,'fro') + norm(A(2:2:end,:) - D,'fro') < 1e-10, ...
    'radgreen interleave');
assert(strcmp(kernel([Ks; Kdp]).sing,'hs'),'radgreen interleave sing');

K = kernel([Ks Kd; Kr Ks]);
A = K.eval(src,targ);
err = norm(A(1:2:end,1:2:end) - S,'fro') + norm(A(1:2:end,2:2:end) - D,'fro') ...
    + norm(A(2:2:end,1:2:end) - Sp,'fro') + norm(A(2:2:end,2:2:end) - S,'fro');
assert(err < 1e-10,'radgreen block interleave');

assert(chk(kernel('radgreen',green,'d'),D),'kernel radgreen dispatch');
assert(chk(kernel('rg',green,'d'),D),'kernel rg dispatch');
assert(chk(kernel('radpp',green,[1e-6 6],'s'),S),'kernel radpp dispatch');
assert(chk(kernel('rp',green,[1e-6 6],'s'),S),'kernel rp dispatch');

end


function test_options()
% name-value parameters, default parameters and the fit options

[src,targ] = geom();
green = @(s,t) chnk.lap2d.green(s,t);
fg = @(s,t) gaussgreen(2,s,t);
cf = [0.7 1.3];
nu = 0.3;

A = kernel.radgreen(green,'c',cf).eval(src,targ);
B = kernel.radgreen(green,'c','coefs',cf).eval(src,targ);
assert(isequal(A,B),'radgreen coefs option');
A = kernel.radgreen(fg,'free_bc',nu).eval(src,targ);
B = kernel.radgreen(fg,'free_bc','nu',nu).eval(src,targ);
assert(isequal(A,B),'radgreen nu option');

w = warning('off','all');
A = kernel.radgreen(green,'c').eval(src,targ);
B = kernel.radgreen(green,'c',[1 1]).eval(src,targ);
A2 = kernel.radgreen(green,'trans').eval(src,targ);
B2 = kernel.radgreen(green,'trans',ones(2)).eval(src,targ);
A3 = kernel.radgreen(green,'c2trans',cf).eval(src,targ);
B3 = kernel.radgreen(green,'c2trans',[cf; cf]).eval(src,targ);
warning(w);
assert(isequal(A,B),'radgreen default coefs');
assert(isequal(A2,B2),'radgreen default trans coefs');
assert(isequal(A3,B3),'radgreen c2trans coefs');

Kt = kernel.radpp(fg,[1e-6 6],'s');
Kl = kernel.radpp(fg,[1e-6 6],'s','opts',struct('eps',1e-6));
K1 = kernel.radpp(fg,6,'s');
assert(Kt.params.rmin == 1e-6 && Kt.params.rmax == 6,'radpp rmin rmax');
assert(K1.params.rmin == 1e-12 && K1.params.rmax == 6,'radpp default rmin');
assert(Kl.params.fit.nbin < Kt.params.fit.nbin,'radpp fit options');
B = kernel.radgreen(fg,'s').eval(src,targ);
assert(norm(K1.eval(src,targ) - B,'fro') < 1e-10*norm(B,'fro'),'radpp default rmin');

end


function test_edge()
% coincident points, in-plane derivatives in 3D, failed singularity detection

green = @(s,t) gaussgreen(2,s,t);
p = [0.3; -0.2];
n1 = [0.6; 0.8];
n2 = [-0.28; 0.96];
s0 = struct('r',p,'n',n1);
t0 = struct('r',p,'n',n2);
err = abs(kernel.radgreen(green,'s').eval(s0,t0) - 1) ...
    + abs(kernel.radgreen(green,'sp').eval(s0,t0)) ...
    + abs(kernel.radgreen(green,'dp').eval(s0,t0) - n1.'*n2);
assert(err < 1e-14,'radgreen coincident points: %.2e',err);

ns = 7;  nt = 9;
s3.r = [randn(2,ns); 0.5*randn(1,ns)];
s3.n = randn(2,ns);  s3.n = s3.n./vecnorm(s3.n);
t3.r = [randn(2,nt); zeros(1,nt)];
t3.n = randn(2,nt);  t3.n = t3.n./vecnorm(t3.n);
dx = t3.r(1,:).' - s3.r(1,:);
dy = t3.r(2,:).' - s3.r(2,:);
dz = t3.r(3,:).' - s3.r(3,:);
G3 = exp(-(dx.^2 + dy.^2 + dz.^2)/2);
dt = dx.*t3.n(1,:).' + dy.*t3.n(2,:).';
ds = dx.*s3.n(1,:) + dy.*s3.n(2,:);
tns = t3.n(1,:).'.*s3.n(1,:) + t3.n(2,:).'.*s3.n(2,:);
ref = struct('s',G3,'d',ds.*G3,'sp',-dt.*G3,'dp',(tns - dt.*ds).*G3);
types = fieldnames(ref);
for i = 1:numel(types)
    A = kernel.radgreen(green,types{i}).eval(s3,t3);
    err = norm(A - ref.(types{i}),'fro')/norm(ref.(types{i}),'fro');
    assert(err < 1e-13,'radgreen in-plane %s: %.2e',types{i},err);
end

w = warning('off','CHNK:radpp:singtype');
K = kernel.radgreen(green,@(D,s,t) NaN*D({},{}),'nderiv',0);
warning(w);
assert(strcmp(K.sing,'hs'),'radgreen: failed detection should give hs');

end


function test_radpp_formats()
% radpp with packed and full derivative arrays, and the fit limits

[src,targ] = geom();
full = @(s,t) gaussgreen(2,s,t);
packed = @(s,t) gausspacked(2,s,t);
for type = {'d','dp','sg','cp','c2trans','clamped_eval'}
    a = targs(type{1});
    A = kernel.radpp(full,[1e-6 6],type{1},a{:}).eval(src,targ);
    B = kernel.radpp(packed,[1e-6 6],type{1},a{:}).eval(src,targ);
    assert(norm(A - B,'fro') < 1e-12*norm(A,'fro'),'radpp packed %s',type{1});
end

K = kernel.radpp(full,[1e-3 6],'s');
s0 = struct('r',[0;0]);
assert(isnan(K.eval(s0,struct('r',[1e-4;0]))),'radpp: no NaN below the fit');
assert(isnan(K.eval(s0,struct('r',[7;0]))),'radpp: no NaN above the fit');
assert(~isnan(K.eval(s0,struct('r',[1;0]))),'radpp: NaN inside the fit');

end


function test_warnings()
% warnings for missing coefficients and failed singularity detection

green = @(s,t) gaussgreen(2,s,t);
lastwarn('');
w = warning('off','all');
kernel.radgreen(green,'c');
msg = lastwarn;
warning(w);
assert(contains(msg,'Missing coefficients'),'radgreen: no missing coefs warning');

lastwarn('','');
w = warning('off','all');
kernel.radgreen(green,@(D,s,t) NaN*D({},{}),'nderiv',0);
[~,id] = lastwarn;
warning(w);
assert(strcmp(id,'CHNK:radpp:singtype'),'radgreen: no singtype warning');

end


function test_chunkermat()
% chunkermat with radgreen kernels against the built-in kernels

cparams = [];
cparams.eps = 1e-10;
pref = [];
pref.k = 16;
chnkr = chunkerfunc(@(t) starfish(t,3,0.3),cparams,pref);

green = @(s,t) chnk.lap2d.green(s,t);
for type = {'s','d','sp'}
    A = chunkermat(chnkr,kernel('lap',type{1}));
    B = chunkermat(chnkr,kernel.radgreen(green,type{1}));
    err = norm(A - B,'fro')/norm(A,'fro');
    assert(err < 1e-9,'chunkermat radgreen %s: %.2e',type{1},err);
end

end


function test_bad_inputs()
% bad inputs are errors

green = @(s,t) gaussgreen(2,s,t);
bad = {@() kernel.radgreen(green,'nonsense'), ...
    @() kernel.radgreen(green,'free_bc'), ...
    @() kernel.radgreen(green,'free_bc',0.3,'dim',3), ...
    @() kernel.radgreen(green,'s',1), ...
    @() kernel.radgreen(green,'s','dim',4), ...
    @() kernel.radgreen(green,'s','nonsense',1), ...
    @() kernel.radgreen(@(s,t) chnk.lap2d.green(s,t),'free_bc',0.3), ...
    @() kernel.radgreen(@notradial,'s'), ...
    @() kernel.radgreen(@wronggrad,'sp'), ...
    @() kernel.radgreen(green,@(D,s,t) D({},{})), ...
    @() kernel.radpp(green)};
for k = 1:numel(bad)
    caught = false;
    try
        bad{k}();
    catch
        caught = true;
    end
    assert(caught,'radgreen: bad input %d not rejected',k);
end

end


function [src,targ] = geom()
% sources and targets on circles, with normals, tangents and data

ns = 7;  nt = 9;
th = 2*pi*rand(1,ns);
src.r = [cos(th); sin(th)];
src.d = [-sin(th); cos(th)];
src.d2 = -src.r;
src.n = src.r;
src.data = randn(2,ns);

th = 2*pi*rand(1,nt);
targ.r = 1.7*[cos(th); sin(th)];
targ.d = 1.7*[-sin(th); cos(th)];
targ.d2 = -targ.r;
targ.n = [cos(th); sin(th)];

end


function varargout = gaussgreen(dim,s,t)
% G = exp(-r^2/2) and its target derivatives as full tensors

nt = size(t,2);  ns = size(s,2);
X = cell(1,dim);
for i = 1:dim
    X{i} = t(i,:).' - s(i,:);
end
G = exp(-sum(cat(3,X{:}).^2,3)/2);

varargout{1} = G;
for k = 1:max(nargout,1)-1
    T = zeros([nt ns repmat(dim,1,k)]);
    for m = 1:dim^k
        sub = cell(1,k);
        [sub{:}] = ind2sub([repmat(dim,1,k) 1],m);
        T(:,:,sub{:}) = hermite(X,[sub{:}]).*G;
    end
    varargout{k+1} = T;
end

end


function p = hermite(X,idx)
% derivative of exp(-x^2/2) over exp(-x^2/2) in the directions idx

x = @(a) X{idx(a)};
d = @(a,b) double(idx(a) == idx(b));
switch numel(idx)
    case 1
        p = -x(1);
    case 2
        p = x(1).*x(2) - d(1,2);
    case 3
        p = -x(1).*x(2).*x(3) + d(1,2)*x(3) + d(1,3)*x(2) + d(2,3)*x(1);
    case 4
        p = x(1).*x(2).*x(3).*x(4) - (d(1,2)*x(3).*x(4) + ...
            d(1,3)*x(2).*x(4) + d(1,4)*x(2).*x(3) + d(2,3)*x(1).*x(4) + ...
            d(2,4)*x(1).*x(3) + d(3,4)*x(1).*x(2)) + ...
            d(1,2)*d(3,4) + d(1,3)*d(2,4) + d(1,4)*d(2,3);
end

end


function varargout = gausspacked(dim,s,t)
% gaussgreen with packed derivatives

out = cell(1,max(nargout,1));
[out{:}] = gaussgreen(dim,s,t);
varargout = out;
for k = 2:nargout-1
    T = out{k+1};
    tup = nchoosek(1:dim+k-1,k) - (0:k-1);
    P = zeros(size(T,1),size(T,2),size(tup,1));
    for m = 1:size(tup,1)
        sub = num2cell(tup(m,:));
        P(:,:,m) = T(:,:,sub{:});
    end
    varargout{k+1} = P;
end

end


function varargout = flexgreen(zk,s,t)
% flexural wave Green's function and its target derivatives

g1 = @(s,t) chnk.flex2d.helmdiffgreen(zk,s,t);
g2 = @(s,t) chnk.flex2d.helmdiffgreen(1i*zk,s,t);
g = chnk.radpp.sumgreens([1 -1]/(2*zk^2),g1,g2);
out = cell(1,max(nargout,1));
[out{:}] = g(s,t);
varargout = out;

end


function varargout = notradial(s,t)
% an anisotropic gaussian

x = t(1,:).' - s(1,:);  y = t(2,:).' - s(2,:);
varargout{1} = exp(-(x.^2 + 2*y.^2)/2);

end


function varargout = wronggrad(s,t)
% gaussian with a gradient that is not that of the value

out = cell(1,max(nargout,1));
[out{:}] = gaussgreen(2,s,t);
varargout = out;
if nargout >= 2
    varargout{2} = 1.1*out{2};
end

end


function v = tcon(T,vs)
% full tensor T (nt,ns,dim,...,dim) contracted with the vectors vs

k = numel(vs);
dim = size(vs{1},3);
v = 0;
for m = 1:dim^k
    sub = cell(1,k);
    [sub{:}] = ind2sub([repmat(dim,1,k) 1],m);
    term = T(:,:,sub{:});
    for j = 1:k
        term = term.*vs{j}(:,:,sub{j});
    end
    v = v + term;
end

end


function a = targs(type)
% the parameter of the types that take one

a = {};
if any(strcmp(type,{'c','cp','cg','trans_rep','trans_rep_prime','trans_rep_grad'}))
    a = {[0.7 1.3i]};
elseif any(strcmp(type,{'c2trans','trans'}))
    a = {[1 2i; 3 4]};
end

end
