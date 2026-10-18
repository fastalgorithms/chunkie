chunkerkernevalmat_chunkgraphTest0();


function chunkerkernevalmat_chunkgraphTest0()
%CHUNKERKERNEVALMAT_CHUNKGRAPHTEST check that chunkerkernevalmat accepts a
% chunkgraph as the target object, and that the result matches passing the
% same target points as a struct.

% source: a circle
cparams = [];
cparams.eps = 1.0e-10;
pref = [];
pref.k = 16;
chnkr = chunkerfunc(@(t) chnk.curves.bymode(t, 3), cparams, pref);

% target: a Y-shaped chunkgraph inside the circle
verts = [0 0; 0 1; -0.7 1.7; 0.7 1.7].';
edge2verts = [1 2 2; 2 3 4];
cgrph = chunkgraph(verts, edge2verts);
cgrph = makedatarows(cgrph,1);
cgrph.data(:) = cgrph.r(1,:).*cgrph.r(2,:);

kern =  (@(t) reshape(t.data(1,:),1,1,[])) * kernel('lap', 's');

targinfo = [];
targinfo.r = cgrph.r(:,:);
targinfo.d = cgrph.d(:,:);
targinfo.d2 = cgrph.d2(:,:);
targinfo.n = cgrph.n(:,:);
targinfo.data = cgrph.data(:,:);

mat_cgrph = chunkerkernevalmat(chnkr, kern, cgrph);
mat_struct = chunkerkernevalmat(chnkr, kern, targinfo);

assert(isequal(size(mat_cgrph), [cgrph.npt, chnkr.npt]));
assert(norm(mat_cgrph - mat_struct, 'fro') == 0);

dens = cos(chnkr.r(1,:)+chnkr.r(2,:)).';
umat = mat_cgrph*dens;
u_cgrph = chunkerkerneval(chnkr,kern,dens,cgrph);
assert(norm(umat - u_cgrph) < 1e-13);

u_struct = chunkerkerneval(chnkr,kern,dens,targinfo);
assert(norm(umat - u_struct) < 1e-13);
end
