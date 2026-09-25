%KERNEL_PARTSTEST check that sums of kernels with different singularities
% keep one part per singularity and that chunkermat uses the matching rule
kernel_partsTest_sum()
kernel_partsTest_scalarops()
kernel_partsTest_mtimes()
kernel_partsTest_interleave()
kernel_partsTest_interleaveops()
kernel_partsTest_chunkermat()
kernel_partsTest_chunkermat_interleave()
kernel_partsTest_chunkermatapply()
kernel_partsTest_rcip()


function kernel_partsTest_sum()
% + and - keep one part per singularity, parts only when types differ

skern = kernel('lap','s');
dkern = kernel('lap','d');
stkern = kernel('lap','stau');
dpkern = kernel('lap','dp');
hskern = kernel('helm','s',1.3);
rkern = kernel('helm1d','s',1.3);
emptykern = kernel(skern.eval);

check_parts(skern + stkern, {'log','pv'}, 'pv');
check_parts(stkern + skern, {'log','pv'}, 'pv');
check_parts(skern - stkern, {'log','pv'}, 'pv');
check_parts(stkern - skern, {'log','pv'}, 'pv');

% same singularity: no parts
check_parts(skern + hskern, {}, 'log');
check_parts(skern - hskern, {}, 'log');

% sums of sums
check_parts((skern + stkern) + (dkern + dpkern), {'hs','log','pv','smooth'}, 'hs');
check_parts((skern + stkern) - (dkern + dpkern), {'hs','log','pv','smooth'}, 'hs');
check_parts((skern + stkern) + hskern, {'log','pv'}, 'pv');
check_parts(hskern - (skern + stkern), {'log','pv'}, 'pv');
check_parts(skern + stkern + dkern + dpkern + hskern, {'hs','log','pv','smooth'}, 'hs');
check_parts((skern + stkern) - (skern + stkern), {'log','pv'}, 'pv');

% removable ranks above smooth
check_parts(dkern + rkern, {'removable','smooth'}, 'removable');
check_parts(rkern + skern, {'log','removable'}, 'log');

% zero kernels do not add a part
check_parts(kernel.zeros() + stkern, {}, 'pv');
check_parts(stkern + kernel.zeros(), {}, 'pv');
check_parts(kernel.zeros() - stkern, {}, 'pv');
check_parts((skern + stkern) + kernel.zeros(), {'log','pv'}, 'pv');

% nan kernels stay nan kernels without parts
knan = kernel.nans() + (skern + stkern);
assert(knan.isnan && isempty(knan.parts));

% a kernel without a singularity type turns off parts
kern = (skern + stkern) + emptykern;
assert(isempty(kern.parts) && strcmp(kern.sing,'pv'));
kern = emptykern + (skern + stkern);
assert(isempty(kern.parts) && isempty(kern.sing));

end


function kernel_partsTest_scalarops()
% scalar operations are applied to every part

skern = kernel('lap','s');
stkern = kernel('lap','stau');
dpkern = kernel('lap','dp');
kern = skern + stkern + dpkern;
keys = {'hs','log','pv'};

check_parts(3*kern, keys, 'hs');
check_parts(kern*3, keys, 'hs');
check_parts(3.*kern, keys, 'hs');
check_parts(kern.*3, keys, 'hs');
check_parts(kern/2, keys, 'hs');
check_parts(kern./2, keys, 'hs');
check_parts(-kern, keys, 'hs');
check_parts(conj(1i*kern), keys, 'hs');
check_parts(2*(-kern)/3 + kern, keys, 'hs');

kzero = 0*kern;
assert(kzero.iszero && isempty(kzero.parts));

end


function kernel_partsTest_mtimes()
% matrix and function handle multiplication is applied to every part

skern = kernel('lap','s');
stkern = kernel('lap','stau');
kern = skern + stkern;
keys = {'log','pv'};

check_parts([2;3]*kern, keys, 'pv');
check_parts(kern*[1 2], keys, 'pv');
check_parts([1 2; 3 4]*([2;3]*kern), keys, 'pv');

hleft = @(t) reshape([1 + t.r(1,:); t.r(2,:)], 2, 1, []);
hright = @(s) reshape([s.r(1,:), 1 + s.r(2,:)], 1, 2, []);
hpoint = @(t) reshape(2 + t.r(1,:), 1, 1, []);
check_parts(hleft*kern, keys, 'pv');
check_parts(kern*hright, keys, 'pv');
check_parts(hpoint*kern, keys, 'pv');
check_parts(kern*hpoint, keys, 'pv');
check_parts(hleft*kern + [1;1]*skern, keys, 'pv');

end


function kernel_partsTest_interleave()
% interleaving keeps one part per singularity

skern = kernel('lap','s');
dkern = kernel('lap','d');
stkern = kernel('lap','stau');
dpkern = kernel('lap','dp');
hskern = kernel('helm','s',1.3);
zkern = kernel.zeros();
emptykern = kernel(skern.eval);

% different types in different blocks
check_parts(kernel([skern, stkern; dkern, 2*skern]), {'log','pv','smooth'}, 'pv');
check_parts(kernel([skern; stkern]), {'log','pv'}, 'pv');
check_parts(kernel([dpkern, skern]), {'hs','log'}, 'hs');

% same type in every block: no parts
check_parts(kernel([skern, hskern; 2*skern, skern]), {}, 'log');

% zero blocks do not add a part
check_parts(kernel([skern, zkern; zkern, stkern]), {'log','pv'}, 'pv');
check_parts(kernel([skern, zkern; zkern, skern]), {}, 'log');

% blocks that already have parts
check_parts(kernel([skern + stkern, dkern; stkern, dpkern]), ...
    {'hs','log','pv','smooth'}, 'hs');
check_parts(kernel([skern + stkern, skern; stkern, skern + stkern]), ...
    {'log','pv'}, 'pv');

% blocks with and without shifted_eval
pvshifted = stkern;
pvshifted.shifted_eval = @(s,t,o) 2*stkern.eval(s,t);
kern = kernel([skern, pvshifted; dkern, skern]);
assert(isa(kern.shifted_eval, 'function_handle'));
check_parts(kern, {'log','pv','smooth'}, 'pv');
check_parts(kernel([skern + pvshifted, dkern]), {'log','pv','smooth'}, 'pv');

% nested interleaving
kint = kernel([skern, stkern; dkern, skern]);
check_parts(kernel([kint, kint]), {'log','pv','smooth'}, 'pv');
check_parts(kernel([kint; kernel([dpkern, zkern])]), {'hs','log','pv','smooth'}, 'hs');

% a block without a singularity type turns off parts
kern = kernel([skern, stkern; emptykern, skern]);
assert(isempty(kern.parts));

end


function kernel_partsTest_interleaveops()
% operations on interleaved kernels

skern = kernel('lap','s');
dkern = kernel('lap','d');
stkern = kernel('lap','stau');
dpkern = kernel('lap','dp');

kint1 = kernel([skern, stkern; dkern, skern]);
kint2 = kernel([stkern, skern; skern, dpkern]);
keys = {'log','pv','smooth'};

check_parts(2*kint1, keys, 'pv');
check_parts(-kint1, keys, 'pv');
check_parts(conj(1i*kint1), keys, 'pv');
check_parts(kint1/3, keys, 'pv');
check_parts([1 2]*kint1, keys, 'pv');
check_parts(kint1*[1; 2], keys, 'pv');
check_parts([1 2; 3 4]*kint1*[0 1; 1 0], keys, 'pv');
check_parts(kint1 + kint2, {'hs','log','pv','smooth'}, 'hs');
check_parts(kint1 - kint2, {'hs','log','pv','smooth'}, 'hs');
check_parts(kint1 + kernel([skern, skern; skern, skern]), keys, 'pv');

end


function kernel_partsTest_chunkermat()
% chunkermat of an operation on K1+K2 equals the sum over K1 and K2 when
% K1 and K2 have different singularities

chnkr = chunkerfunc(@(t) starfish(t,3,0.1));

skern = kernel('lap','s');
dkern = kernel('lap','d');
stkern = kernel('lap','stau');
dpkern = kernel('lap','dp');
hleft = @(t) reshape([1 + t.r(1,:); t.r(2,:)], 2, 1, []);

ops = {@(k) k, @(k) -k, @(k) conj(1i*k), @(k) 3*k, @(k) k/2, ...
    @(k) [2;3]*k, @(k) k*[1 2], @(k) hleft*k};
pairs = {{skern, stkern}, {skern + dkern, stkern}, {skern, dpkern}, ...
    {dkern, stkern}};
optslist = {struct(), struct('nonsmoothonly',true), struct('corrections',true)};

for iop = 1:numel(ops)
    for ipair = 1:numel(pairs)
        k1 = ops{iop}(pairs{ipair}{1});
        k2 = ops{iop}(pairs{ipair}{2});
        kern = ops{iop}(pairs{ipair}{1} + pairs{ipair}{2});
        for iopts = 1:numel(optslist)
            A = chunkermat(chnkr, kern, optslist{iopts});
            B = chunkermat(chnkr, k1, optslist{iopts}) + ...
                chunkermat(chnkr, k2, optslist{iopts});
            assert(norm(full(A-B),'fro') < 1e-12*norm(full(B),'fro'));
        end
    end
end

% skipped panels, their self blocks are left as NaN
ilist = [1; 2; 5];
kern = skern + stkern + dpkern;
A = chunkermat(chnkr, kern, [], ilist);
B = chunkermat(chnkr, skern, [], ilist) + chunkermat(chnkr, stkern, [], ilist) ...
    + chunkermat(chnkr, dpkern, [], ilist);
assert(isequal(isnan(A), isnan(B)));
ifinite = ~isnan(B);
assert(norm(A(ifinite)-B(ifinite)) < 1e-12*norm(B(ifinite)));

% reusing the quadrature rules stored in opts
[A, opts] = chunkermat(chnkr, kern);
assert(all(isfield(opts.auxquads, {'ggqlog','ggqpv','ggqhs'})));
B = chunkermat(chnkr, kern, opts);
assert(norm(A-B,'fro') == 0);

% stored rules for a different order are not reused
chnkr8 = chunkerfunc(@(t) starfish(t), struct('maxchunklen',1), struct('k',8));
assert(chnkr8.k == 8 && chnkr.k ~= 8);
A = chunkermat(chnkr8, kern);
B = chunkermat(chnkr8, kern, opts);
assert(norm(A-B,'fro') == 0);

end


function kernel_partsTest_chunkermat_interleave()
% chunkermat of an interleaved kernel equals the blockwise matrices

chnkr = chunkerfunc(@(t) starfish(t,3,0.1));
n = chnkr.npt;

skern = kernel('lap','s');
dkern = kernel('lap','d');
stkern = kernel('lap','stau');
dpkern = kernel('lap','dp');

blockslist = {[skern, stkern; dkern, 2*skern], ...
    [skern + stkern, dkern; stkern, dpkern], ...
    [skern, kernel.zeros(); kernel.zeros(), stkern]};

for iblocks = 1:numel(blockslist)
    blocks = blockslist{iblocks};
    B = zeros(2*n);
    for i = 1:2
        for j = 1:2
            B(i:2:end,j:2:end) = chunkermat(chnkr, blocks(i,j));
        end
    end

    A = chunkermat(chnkr, kernel(blocks));
    assert(norm(A-B,'fro') < 1e-12*norm(B,'fro'));

    A = chunkermat(chnkr, 2*kernel(blocks));
    assert(norm(A-2*B,'fro') < 1e-12*norm(B,'fro'));

    A = chunkermat(chnkr, [1 2; 3 4]*kernel(blocks));
    C = kron(speye(n), [1 2; 3 4])*B;
    assert(norm(A-C,'fro') < 1e-12*norm(C,'fro'));
end

end


function kernel_partsTest_chunkermatapply()
% chunkermatapply uses the parts through chunkermat

chnkr = chunkerfunc(@(t) starfish(t,3,0.1));

kern = kernel('lap','s') + kernel('lap','dp');
dens = randn(chnkr.npt,1);

A = chunkermat(chnkr, kern);
u = chunkermatapply(chnkr, kern, dens);
assert(norm(u - A*dens) < 1e-10*norm(A*dens));

end


function kernel_partsTest_rcip()
% RCIP uses the shifted evaluator of each part

verts = [0 1 1; 0 0 1];
cg = chunkgraph(verts, [1:3; [2:3 1]]);
nedge = size(cg.echnks,2);

skern = kernel('lap','s');
stkern = kernel('lap','stau');

% position dependent pv kernel, so a missing shift gives the wrong answer
wfun = @(t,o) (1 + sum((t.r(:,:) + o(:)).^2, 1)).';
pvkern = stkern;
pvkern.eval = @(s,t) wfun(t,[0;0]).*stkern.eval(s,t);
pvkern.shifted_eval = @(s,t,o) wfun(t,o).*stkern.eval(s,t);

kern = skern + pvkern;
kernsingle = kern;
kernsingle.parts = [];
B = chunkermat(cg, kernsingle);

A = chunkermat(cg, kern);
assert(norm(A-B,'fro') < 1e-8*norm(B,'fro'));

% one kernel per pair of edges
kernmat(nedge,nedge) = kernel();
for i = 1:nedge
    for j = 1:nedge
        kernmat(i,j) = kern;
    end
end
A = chunkermat(cg, kernmat);
assert(norm(A-B,'fro') < 1e-8*norm(B,'fro'));

% interleaved kernel where only some blocks have shifted_eval
kint = kernel([skern, pvkern; pvkern, skern]);
check_parts(kint, {'log','pv'}, 'pv');
kintsingle = kint;
kintsingle.parts = [];
A = chunkermat(cg, kint);
B = chunkermat(cg, kintsingle);
assert(norm(A-B,'fro') < 1e-8*norm(B,'fro'));

end


function check_parts(kern, keys, sing)
% check the parts of kern against the expected keys and singularity

assert(strcmp(kern.sing, sing));
if isempty(keys)
    assert(isempty(kern.parts));
    return
end
assert(isequal(sort(fieldnames(kern.parts)), sort(keys(:))));

src = []; src.r = randn(2,5); src.d = randn(2,5); src.d2 = randn(2,5);
src.n = randn(2,5);
targ = []; targ.r = randn(2,4); targ.d = randn(2,4); targ.d2 = randn(2,4);
targ.n = randn(2,4);
o = randn(2,1);

val = kern.eval(src,targ);
partsum = 0;
shiftedpartsum = 0;
for sname = fieldnames(kern.parts)'
    part = kern.parts.(sname{1});
    assert(isempty(part.parts));
    assert(strcmp(part.sing, sname{1}));
    assert(isequal(part.opdims, kern.opdims));
    partsum = partsum + part.eval(src,targ);
    if isa(part.shifted_eval, 'function_handle')
        shiftedpartsum = shiftedpartsum + part.shifted_eval(src,targ,o);
    else
        shiftedpartsum = shiftedpartsum + part.eval(src,targ);
    end
end
assert(norm(partsum - val) <= 1e-12*max(1, norm(val)));

if isa(kern.shifted_eval, 'function_handle')
    val = kern.shifted_eval(src,targ,o);
    assert(norm(shiftedpartsum - val) <= 1e-12*max(1, norm(val)));
end

end
