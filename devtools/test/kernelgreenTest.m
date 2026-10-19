%KERNELGREENTEST test kernel.green and kernel.radpp

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
test_nonsymmetric();
test_fmm();
test_grad();


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

% noise in f below the tolerance: a warning, not maxsub
f = @(r) exp(1i*100*r)./r;
lastwarn('');
warning('off','CHNK:radpp:fit:noise');
info = chnk.radpp.fit(f,1e-3,10,struct('eps',1e-15));
warning('on','CHNK:radpp:fit:noise');
[~,id] = lastwarn();
assert(strcmp(id,'CHNK:radpp:fit:noise'),'radpp fit: no noise warning');
assert(info.err < 1e-11,'radpp fit: noise level %.1e',info.err);

% a sharp f stalls below 1e-4 before it is resolved, unless stopnoise is
% false
f = @(r) tanh(100*(r-1));
r = linspace(1e-3,10,10000).';
lastwarn('');
warning('off','CHNK:radpp:fit:noise');
info = chnk.radpp.fit(f,1e-3,10,struct('eps',1e-13));
warning('on','CHNK:radpp:fit:noise');
[~,id] = lastwarn();
assert(strcmp(id,'CHNK:radpp:fit:noise'),'radpp fit: sharp f should stall');
lastwarn('');
info = chnk.radpp.fit(f,1e-3,10,struct('eps',1e-13,'stopnoise',false));
assert(isempty(lastwarn()),'radpp fit: stopnoise = false warned');
err = max(abs(chnk.radpp.fiteval(r,info.breaks,info.coefs) - f(r)));
assert(err < 1e-12,'radpp fit: stopnoise = false, error %.1e',err);

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
        C = [C, chnk.green.unpack(out{k+1},dim,k)];
    end

    pout = cell(1,5);
    [pout{:}] = chnk.radpp.green(@(r) gaussprof(r),dim,zeros(dim,1),X(:,1:dim).');
    P = pout{1}(:);
    for k = 1:4
        P = [P, reshape(pout{k+1},n,[])];
    end
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
    A = kernel.green(full,types{i},a{:}).eval(src,targ);
    B = kernel.green(@(s,t) gaussgreen(2,s,t,true),types{i},a{:}).eval(src,targ);
    err = norm(A - B,'fro')/norm(A,'fro');
    assert(err < 1e-14,'green packed format %s: %.2e',types{i},err);
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
    K = kernel.green(green,tab{i,1},a{:});
    A = K.eval(src,targ);
    B = chnk.lap2d.kern(src,targ,tab{i,2},a{:});
    err = norm(A - B,'fro')/norm(B,'fro');
    assert(isequal(size(A),size(B)),'green lap %s: wrong size',tab{i,1});
    assert(err < 1e-11,'green lap %s: %.2e',tab{i,1},err);
end

% the two row kernels are the rows of s, sp and d, dp
S = kernel.green(green,'s').eval(src,targ);
Sp = kernel.green(green,'sp').eval(src,targ);
D = kernel.green(green,'d').eval(src,targ);
Dp = kernel.green(green,'dp').eval(src,targ);
A = kernel.green(green,'s2trans').eval(src,targ);
assert(norm(A(1:2:end,:) - S,'fro') + norm(A(2:2:end,:) - Sp,'fro') < 1e-12, ...
    'green s2trans');
A = kernel.green(green,'d2trans').eval(src,targ);
assert(norm(A(1:2:end,:) - D,'fro') + norm(A(2:2:end,:) - Dp,'fro') < 1e-12, ...
    'green d2trans');


% custom kernels
cfun = @(D,s,t) cf(1)*D({t.n},{s.n}) + cf(2)*D({t.n},{});
A = kernel.green(green,cfun,'nderiv',2).eval(src,targ);
B = kernel.green(green,'cp',cf).eval(src,targ);
assert(norm(A - B,'fro') < 1e-12*norm(B,'fro'),'green custom scalar');
sfun = @(D,s,t) cat(3,D({},{}),D({t.n},{}));
K = kernel.green(green,sfun,'nderiv',1,'opdims',[2 1]);
A = K.eval(src,targ);
B = kernel.green(green,'s2trans').eval(src,targ);
assert(norm(A - B,'fro') < 1e-12*norm(B,'fro'),'green custom vector');
assert(strcmp(K.sing,'log'),'green custom sing');

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
    K = kernel.green(green,tab{i,1},a{:});
    A = K.eval(src,targ);
    B = chnk.helm2d.kern(zk,src,targ,tab{i,2},a{:});
    err = norm(A - B,'fro')/norm(B,'fro');
    assert(isequal(size(A),size(B)),'green helm %s: wrong size',tab{i,1});
    assert(err < 1e-11,'green helm %s: %.2e',tab{i,1},err);
    assert(isequal(K.opdims,[size(A,1)/numel(targ.n(1,:)), ...
        size(A,2)/numel(src.n(1,:))]),'green helm %s: opdims',tab{i,1});
end

end


function test_flex()
% plate kernels against chnk.flex2d.kern

[src,targ] = geom();
zk = 1.3;
nu = 0.3;
green = chnk.green.sum([1 -1]/(2*zk^2), ...
    @(s,t) chnk.flex2d.helmdiffgreen(zk,s,t), ...
    @(s,t) chnk.flex2d.helmdiffgreen(1i*zk,s,t));

tab = {'clamped_bc','clamped_plate_bcs',0; ...
    'clamped_eval','clamped_plate_eval',0; ...
    'supported_bc','supported_plate_bcs',1; ...
    'supported_eval','supported_plate_eval',1; ...
    'free_bc','free_plate_bcs',1; 'free_eval','free_plate_eval',1};
for i = 1:size(tab,1)
    a = {};
    if tab{i,3}, a = {nu}; end
    A = kernel.green(green,tab{i,1},a{:}).eval(src,targ);
    B = chnk.flex2d.kern(zk,src,targ,tab{i,2},a{:});
    err = norm(A - B,'fro')/norm(B,'fro');
    assert(isequal(size(A),size(B)),'green flex %s: wrong size',tab{i,1});
    assert(err < 1e-9,'green flex %s: %.2e',tab{i,1},err);
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
    K = kernel.green(green,types{i},a{:},'dim',3);
    A = K.eval(src,targ);
    err = norm(A - ref.(types{i}),'fro')/norm(ref.(types{i}),'fro');
    assert(isequal(size(A),size(ref.(types{i}))),'green 3d %s: size',types{i});
    assert(err < 1e-12,'green 3d %s: %.2e',types{i},err);
end

% packed derivative arrays
packed = @(s,t) gaussgreen(3,s,t,true);
for type = {'dp','cg','trans'}
    a = targs(type{1});
    A = kernel.green(green,type{1},a{:},'dim',3).eval(src,targ);
    B = kernel.green(packed,type{1},a{:},'dim',3).eval(src,targ);
    assert(norm(A - B,'fro') < 1e-14*norm(A,'fro'),'green 3d packed %s',type{1});
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
        Kg = kernel.green(greens{j},types{i},a{:});
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
    Kg = kernel.green(green,type{1},a{:},'dim',3);
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
gs = chnk.green.sum([2 -1],g1,g2);
A = kernel.green(gs,'sp').eval(src,targ);
B = 2*kernel.green(g1,'sp').eval(src,targ) - ...
    kernel.green(g2,'sp').eval(src,targ);
assert(norm(A - B,'fro') < 1e-12*norm(B,'fro'),'green: sumgreens');

end


function test_sing()
% singularity types of the kernels

green = @(s,t) chnk.lap2d.green(s,t);
tab = {'s','log'; 'sp','log'; 'd','log'; 'dp','hs'; 'sg','pv'; 'dg','hs'; ...
    'st','pv'; 'hilb','pv'};
for i = 1:size(tab,1)
    K = kernel.green(green,tab{i,1});
    assert(strcmp(K.sing,tab{i,2}),'green sing %s: %s',tab{i,1},K.sing);
end

green = @(s,t) chnk.helm2d.green(1.1,s,t);
tab = {'s','log'; 'sp','log'; 'd','log'; 'dp','hs'; 'c','log'; 'cp','hs'; ...
    'sg','pv'; 'dg','hs'; 'cg','hs'; 'c2trans','hs'; 'trans','hs'; ...
    'trans_rep','log'; 'trans_rep_prime','hs'; 'trans_rep_grad','hs'};
for i = 1:size(tab,1)
    a = targs(tab{i,1});
    K = kernel.green(green,tab{i,1},a{:});
    assert(strcmp(K.sing,tab{i,2}),'green helm sing %s: %s',tab{i,1},K.sing);
end

K = kernel.green(green,'s','sing','hs');
assert(strcmp(K.sing,'hs'),'green: sing option');

zk = 1.3;
green = chnk.green.sum([1 -1]/(2*zk^2), ...
    @(s,t) chnk.flex2d.helmdiffgreen(zk,s,t), ...
    @(s,t) chnk.flex2d.helmdiffgreen(1i*zk,s,t));
K = kernel.green(green,'clamped_bc');
assert(strcmp(K.sing,'log'),'green sing clamped_bc: %s',K.sing);

end


function test_arithmetic()
% kernel arithmetic, interleaving and dispatch with green and radpp

[src,targ] = geom();
green = @(s,t) chnk.lap2d.green(s,t);
Ks = kernel.green(green,'s');
Kd = kernel.green(green,'d');
Kdp = kernel.green(green,'dp');
Kr = kernel.radpp(green,[1e-6 6],'sp');
Kb = kernel('lap','d');
S = Ks.eval(src,targ);
D = Kd.eval(src,targ);
Sp = Kr.eval(src,targ);
a = 2 + 1i;

chk = @(K,B) norm(K.eval(src,targ) - B,'fro') < 1e-10*norm(B,'fro');
assert(chk(Ks + Kd,S + D),'green plus');
assert(chk(Ks - Kd,S - D),'green minus');
assert(chk(-Ks,-S),'green uminus');
assert(chk(a.*Ks,a*S),'green times');
assert(chk(Ks/a,S/a),'green divide');
assert(chk(Ks + Kb,S + Kb.eval(src,targ)),'green plus built-in');
assert(chk(Ks + Kr,S + Sp),'green plus radpp');
Ksum = Ks + Kd;
assert(strcmp(Ksum.sing,'log'),'green plus sing');
Ksum = Ks + Kdp;
assert(strcmp(Ksum.sing,'hs'),'green plus sing hs');

K = kernel([Ks; Kd]);
A = K.eval(src,targ);
assert(isequal(K.opdims,[2 1]),'green interleave opdims');
assert(norm(A(1:2:end,:) - S,'fro') + norm(A(2:2:end,:) - D,'fro') < 1e-10, ...
    'green interleave');
assert(strcmp(kernel([Ks; Kdp]).sing,'hs'),'green interleave sing');

K = kernel([Ks Kd; Kr Ks]);
A = K.eval(src,targ);
err = norm(A(1:2:end,1:2:end) - S,'fro') + norm(A(1:2:end,2:2:end) - D,'fro') ...
    + norm(A(2:2:end,1:2:end) - Sp,'fro') + norm(A(2:2:end,2:2:end) - S,'fro');
assert(err < 1e-10,'green block interleave');

assert(chk(kernel('green',green,'d'),D),'kernel green dispatch');

% sums and multiples of radpp kernels of one type are radpp kernels
g1 = @(s,t) chnk.helm2d.green(1.1,s,t);
g2 = @(s,t) chnk.helm2d.green(2.3,s,t);
K1 = kernel.radpp(g1,[1e-6 6],'c',[1 2i]);
K2 = kernel.radpp(g2,[1e-5 8],'c',[1 2i]);
A1 = K1.eval(src,targ);  A2 = K2.eval(src,targ);
tab = {2*K1 - K2, 2*A1 - A2; K1 + K2/3, A1 + A2/3; -K1, -A1; K1.*1i, 1i*A1};
for i = 1:size(tab,1)
    K = tab{i,1};
    assert(strcmp(K.name,'radpp'),'radpp sum %d: not radpp',i);
    assert(norm(K.eval(src,targ) - tab{i,2},'fro') < 1e-13*norm(tab{i,2},'fro'), ...
        'radpp sum %d',i);
end
K = K1 + K2;
assert(K.params.rmin == 1e-5 && K.params.rmax == 6,'radpp sum: range');
K = K1 + kernel.radpp(g1,6,'sp');
assert(~strcmp(K.name,'radpp'),'radpp sum: different types merged');

% the fmms of sums, scalings and interleaved kernels, with two densities
if exist(['fmm2d.' mexext],'file') == 3
    zk = 1.1;
    gh = @(s,t) chnk.helm2d.green(zk,s,t);
    fmm = @(eps,s,pg,t,pgt) hfmm2d(eps,zk,s,pg,t,pgt);
    Ks = kernel.green(gh,'s','fmm',fmm);
    Kd = kernel.green(gh,'d','fmm',fmm);
    Kr = kernel.radpp(gh,[1e-6 6],'sp','fmm',fmm);
    Kg = kernel.green(chnk.green.sum([2 -1],gh,gh),'s','fmm', ...
        @(eps,s,pg,t,pgt) fmm(eps,s,pg,t,pgt));
    kerns = {Ks + Kd, Ks - Kr, -Kd, a.*Ks, Ks./a, Ks/a, conj(Kd), ...
        kernel([Ks Kd; Kr Ks]), kernel([Ks; Kr]), Kg};
    for i = 1:numel(kerns)
        K = kerns{i};
        A = K.eval(src,targ);
        sigma = randn(size(A,2),2);
        u = K.fmm(1e-12,src,targ,sigma);
        assert(norm(u - A*sigma,'fro') < 1e-9*norm(A*sigma,'fro'), ...
            'green fmm arithmetic %d',i);
    end
end
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

A = kernel.green(green,'c',cf).eval(src,targ);
B = kernel.green(green,'c','coefs',cf).eval(src,targ);
assert(isequal(A,B),'green coefs option');
A = kernel.green(fg,'free_bc',nu).eval(src,targ);
B = kernel.green(fg,'free_bc','nu',nu).eval(src,targ);
assert(isequal(A,B),'green nu option');

w = warning('off','all');
A = kernel.green(green,'c').eval(src,targ);
B = kernel.green(green,'c',[1 1]).eval(src,targ);
A2 = kernel.green(green,'trans').eval(src,targ);
B2 = kernel.green(green,'trans',ones(2)).eval(src,targ);
A3 = kernel.green(green,'c2trans',cf).eval(src,targ);
B3 = kernel.green(green,'c2trans',[cf; cf]).eval(src,targ);
warning(w);
assert(isequal(A,B),'green default coefs');
assert(isequal(A2,B2),'green default trans coefs');
assert(isequal(A3,B3),'green c2trans coefs');

Kt = kernel.radpp(fg,[1e-6 6],'s');
Kl = kernel.radpp(fg,[1e-6 6],'s','opts',struct('eps',1e-6));
K1 = kernel.radpp(fg,6,'s');
assert(Kt.params.rmin == 1e-6 && Kt.params.rmax == 6,'radpp rmin rmax');
assert(K1.params.rmin == 1e-12 && K1.params.rmax == 6,'radpp default rmin');
assert(Kl.params.nbin < Kt.params.nbin,'radpp fit options');
B = kernel.green(fg,'s').eval(src,targ);
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
err = abs(kernel.green(green,'s').eval(s0,t0) - 1) ...
    + abs(kernel.green(green,'sp').eval(s0,t0)) ...
    + abs(kernel.green(green,'dp').eval(s0,t0) - n1.'*n2);
assert(err < 1e-14,'green coincident points: %.2e',err);

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
    A = kernel.radpp(green,[1e-6 8],types{i}).eval(s3,t3);
    err = norm(A - ref.(types{i}),'fro')/norm(ref.(types{i}),'fro');
    assert(err < 1e-10,'radpp in-plane %s: %.2e',types{i},err);
end

w = warning('off','CHUNKIE:kernel:green:sing');
K = kernel.green(green,@(D,s,t) NaN*D({},{}),'nderiv',0);
warning(w);
assert(strcmp(K.sing,'hs'),'green: failed detection should give hs');

end


function test_radpp_formats()
% radpp with packed and full derivative arrays, and the fit limits

[src,targ] = geom();
full = @(s,t) gaussgreen(2,s,t);
packed = @(s,t) gaussgreen(2,s,t,true);
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
kernel.green(green,'c');
msg = lastwarn;
warning(w);
assert(contains(msg,'Missing coefficients'),'green: no missing coefs warning');

lastwarn('','');
w = warning('off','all');
kernel.green(green,@(D,s,t) NaN*D({},{}),'nderiv',0);
[~,id] = lastwarn;
warning(w);
assert(strcmp(id,'CHUNKIE:kernel:green:sing'),'green: no singtype warning');

end


function test_chunkermat()
% chunkermat with green kernels against the built-in kernels

cparams = [];
cparams.eps = 1e-6;
pref = [];
pref.k = 16;
chnkr = chunkerfunc(@(t) starfish(t,3,0.3),cparams,pref);

green = @(s,t) chnk.lap2d.green(s,t);
for type = {'s','d','sp'}
    A = chunkermat(chnkr,kernel('lap',type{1}));
    B = chunkermat(chnkr,kernel.green(green,type{1}));
    err = norm(A - B,'fro')/norm(A,'fro');
    assert(err < 1e-9,'chunkermat green %s: %.2e',type{1},err);
end

end


function test_bad_inputs()
% bad inputs are errors

green = @(s,t) gaussgreen(2,s,t);
bad = {@() kernel.green(green,'nonsense'), ...
    @() kernel.green(green,'free_bc'), ...
    @() kernel.green(green,'free_bc',0.3,'dim',3), ...
    @() kernel.green(green,'s',1), ...
    @() kernel.green(green,'s','dim',4), ...
    @() kernel.green(green,'s','nonsense',1), ...
    @() kernel.green(@(s,t) chnk.lap2d.green(s,t),'free_bc',0.3), ...
    @() kernel.radpp(@(s,t) exp(-((t(1,:).'-s(1,:)).^2 + ...
        2*(t(2,:).'-s(2,:)).^2)/2),[1e-6 6],'s'), ...
    @() kernel.green(@wronggrad,'sp'), ...
    @() kernel.green(green,@(D,s,t) D({},{})), ...
    @() kernel.radpp(green)};
for k = 1:numel(bad)
    caught = false;
    try
        bad{k}();
    catch
        caught = true;
    end
    assert(caught,'green: bad input %d not rejected',k);
end

end


function test_nonsymmetric()
% a Green's function that is not symmetric: G(x,y) = g(x-y) with g an
% anisotropic, shifted gaussian

[src,targ] = geom();
A = [1.3 0.4; 0.4 0.7];  c = [0.3; -0.2];
green = @(s,t) shiftgauss(A,c,s,t);
K = kernel.green(green,'dp');
X = permute(targ.r,[2 3 1]) - permute(src.r,[3 2 1]) - permute(c,[3 2 1]);
q = sum(X.*sum(permute(A,[3 4 1 2]).*permute(X,[1 2 4 3]),4),3);
g = exp(-q/2);
AX = sum(permute(A,[3 4 1 2]).*permute(X,[1 2 4 3]),4);
% D' = -n_x' H n_y with H the hessian of g
tn = permute(targ.n,[2 3 1]);  sn = permute(src.n,[3 2 1]);
H = @(i,j) (AX(:,:,i).*AX(:,:,j) - A(i,j)).*g;
ref = 0;
for i = 1:2
    for j = 1:2
        ref = ref - tn(:,:,i).*H(i,j).*sn(:,:,j);
    end
end
B = K.eval(src,targ);
assert(norm(B - ref,'fro') < 1e-13*norm(ref,'fro'),'green nonsymmetric dp');

caught = false;
try
    kernel.radpp(green,[1e-6 6],'s');
catch
    caught = true;
end
assert(caught,'radpp: nonradial green not rejected');

end


function test_fmm()
% the generic fmm against direct evaluation, with several densities

if exist(['fmm2d.' mexext],'file') ~= 3
    return
end
[src,targ] = geom();
ns = 150;  nt = 120;  ndens = 2;
th = 2*pi*rand(1,ns);
src.r = [cos(th); sin(th)];  src.n = src.r;
th = 2*pi*rand(1,nt);
targ.r = 1.7*[cos(th); sin(th)];  targ.n = [cos(th); sin(th)];

zk = 1.1;
green = @(s,t) chnk.helm2d.green(zk,s,t);
fmm = @(eps,s,pg,t,pgt) hfmm2d(eps,zk,s,pg,t,pgt);
types = {'s','d','sp','dp','c','cp','sg','dg','cg','st','s2trans', ...
    'd2trans','c2trans','trans','trans_rep','trans_rep_prime', ...
    'trans_rep_grad','clamped_bc'};
for i = 1:numel(types)
    a = targs(types{i});
    K = kernel.green(green,types{i},a{:},'fmm',fmm);
    assert(~isempty(K.fmm),'green fmm %s: missing',types{i});
    A = K.eval(src,targ);
    sigma = randn(size(A,2),ndens) + 1i*randn(size(A,2),ndens);
    u = K.fmm(1e-12,src,targ,sigma);
    err = norm(u - A*sigma,'fro')/norm(A*sigma,'fro');
    assert(err < 1e-9,'green fmm %s: %.2e',types{i},err);
    u1 = K.fmm(1e-12,src,targ,sigma(:,1));
    err = norm(u1 - A*sigma(:,1))/norm(A*sigma(:,1));
    assert(err < 1e-9,'green fmm %s single density: %.2e',types{i},err);
end

K = kernel.green(green,'free_eval',0.3,'fmm',fmm);
assert(isempty(K.fmm),'green fmm: free_eval has no fmm');

% radpp with the fmm of the exact Green's function
K = kernel.radpp(green,[1e-6 6],'c',[1 1i],'fmm',fmm);
A = K.eval(src,targ);
sigma = randn(size(A,2),1);
err = norm(K.fmm(1e-12,src,targ,sigma) - A*sigma)/norm(A*sigma);
assert(err < 1e-9,'radpp fmm: %.2e',err);

end


function test_grad()
% grad and normaltrace (ntrace) against the gradient and normal derivative
% types, and their fmms

[src,targ] = geom();
zk = 1.3;
green = @(s,t) chnk.helm2d.green(zk,s,t);
fmm = @(eps,s,pg,t,pgt) hfmm2d(eps,zk,s,pg,t,pgt);
cf = [0.7 1.3i];
tab = {kernel.green(green,'s','fmm',fmm), 'sg', 'sp'; ...
    kernel.green(green,'d','fmm',fmm), 'dg', 'dp'; ...
    kernel.green(green,'c',cf,'fmm',fmm), 'cg', 'cp'; ...
    kernel.radpp(green,[1e-6 6],'s'), 'sg', 'sp'; ...
    kernel.radpp(green,[1e-6 6],'d'), 'dg', 'dp'; ...
    kernel.green(green,'trans_rep',cf,'fmm',fmm), 'trans_rep_grad', ...
        'trans_rep_prime'; ...
    kernel.radpp(green,[1e-6 6],'trans_rep',cf), 'trans_rep_grad', ...
        'trans_rep_prime'};
for i = 1:size(tab,1)
    K = tab{i,1};
    a = targs(tab{i,2});
    Ag = kernel.green(green,tab{i,2},a{:}).eval(src,targ);
    An = kernel.green(green,tab{i,3},a{:}).eval(src,targ);
    Kg = grad(K);
    Kn = normaltrace(K);
    assert(isequal(Kg.opdims,[2 1].*K.opdims),'grad %d: opdims',i);
    assert(norm(Kg.eval(src,targ) - Ag,'fro') < 1e-10*norm(Ag,'fro'),'grad %d',i);
    assert(norm(Kn.eval(src,targ) - An,'fro') < 1e-10*norm(An,'fro'), ...
        'normaltrace %d',i);
    assert(isequal(ntrace(K).eval(src,targ),Kn.eval(src,targ)),'ntrace %d',i);
    assert(strcmp(Kn.sing,kernel.green(green,tab{i,3},a{:}).sing), ...
        'normaltrace sing %d',i);
    if ~isempty(K.fmm)
        sig = randn(size(Ag,2),2);
        assert(norm(Kg.fmm(1e-12,src,targ,sig) - Ag*sig,'fro') < ...
            1e-10*norm(Ag*sig,'fro'),'grad fmm %d',i);
        assert(norm(Kn.fmm(1e-12,src,targ,sig) - An*sig,'fro') < ...
            1e-10*norm(An*sig,'fro'),'normaltrace fmm %d',i);
    end
end

% the built-in kernels
tab = {kernel.lap2d('s'), 'sp', 'sg'; kernel.lap2d('c',cf), 'cp', 'cg'; ...
    kernel.helm2d('d',zk), 'dp', ''; kernel.helm2d('trans_rep',zk,cf), ...
    'trans_rep_prime', ''};
for i = 1:size(tab,1)
    K = tab{i,1};
    c = {};
    if isfield(K.params,'coefs'), c = {K.params.coefs}; end
    if strcmp(K.name,'laplace')
        mk = @(type) kernel.lap2d(type,c{:});
    else
        mk = @(type) kernel.helm2d(type,zk,c{:});
    end
    A = mk(tab{i,2}).eval(src,targ);
    assert(norm(normaltrace(K).eval(src,targ) - A,'fro') < ...
        1e-14*norm(A,'fro'),'built-in normaltrace %d',i);
    if ~isempty(tab{i,3})
        A = mk(tab{i,3}).eval(src,targ);
        assert(norm(grad(K).eval(src,targ) - A,'fro') < ...
            1e-14*norm(A,'fro'),'built-in grad %d',i);
    end
end

caught = false;
try
    grad(kernel.lap2d('sp'));
catch
    caught = true;
end
assert(caught,'grad: prime type not rejected');

end


function varargout = shiftgauss(A,c,s,t)
% exp(-(x-y-c)'A(x-y-c)/2) and its target derivatives, packed
X = permute(t,[2 3 1]) - permute(s,[3 2 1]) - permute(c,[3 2 1]);
AX = sum(permute(A,[3 4 1 2]).*permute(X,[1 2 4 3]),4);
g = exp(-sum(X.*AX,3)/2);
varargout{1} = g;
varargout{2} = -AX.*g;
varargout{3} = cat(3,(AX(:,:,1).^2 - A(1,1)).*g, ...
    (AX(:,:,1).*AX(:,:,2) - A(1,2)).*g,(AX(:,:,2).^2 - A(2,2)).*g);
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


function varargout = gaussgreen(dim,s,t,packed)
% G = exp(-r^2/2) and its target derivatives, as full tensors or packed

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
    if nargin > 3 && packed && k > 1
        tup = num2cell(nchoosek(1:dim+k-1,k) - (0:k-1),1);
        T = reshape(T,nt*ns,[]);
        T = reshape(T(:,sub2ind([repmat(dim,1,k) 1],tup{:})),nt,ns,[]);
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


function f = gaussprof(r)
% on-axis profiles of exp(-r^2/2), see chnk.radpp.green
r = r(:);
g = exp(-r.^2/2);
f = [g, -r.*g, (r.^2-1).*g, -g, (3*r-r.^3).*g, r.*g, ...
    (r.^4-6*r.^2+3).*g, (1-r.^2).*g, 3*g];
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
