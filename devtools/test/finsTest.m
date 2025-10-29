
verts = [0 0.3 1 1 .5;
         0 1 1 0 0];
evs = [[1;2],[2;3],[3;4],[4;5],[5;1]];
evs = flipud(evs);
cgrph = chunkgraph(verts,evs);
figure(1);clf
plot(cgrph)
quiver(cgrph)

x = linspace(0,1,80)+1e-2; x = x(1:end-1); y = x;
[xx,yy] = meshgrid(x,y);
targs = []; targs.r = [xx(:)';yy(:)'];
dkern = kernel('lap','d');
skern = kernel('lap','s');
spkern = kernel('lap','sp');
dpkern = kernel('lap','dp');

nedge = size(evs,2);
matkerns(nedge,nedge) = kernel();

matkerns(1:3,1:3) = 2*spkern;


matkerns(4:5,1:3) = -2*skern;
matkerns(1:3,4:5) = 2*dpkern;
matkerns(4:5,4:5) = -2*dkern;

% matkerns(1:5,1:5) = skern;
chnkg_n = slicegraph(cgrph,1:3);
chnkg_d = slicegraph(cgrph,4:5);

A = chunkermat(cgrph,matkerns);
A = A + eye(size(A));
source = [1.5, 1.5]';
str = []; str.r = source;
rhs_n = 2*spkern.eval(str,chnkg_n);
rhs_d = -2*skern.eval(str,chnkg_d);
rhs = -[rhs_n; rhs_d];
b = A\rhs;

plotkerns(1,nedge) = kernel();
plotkerns(1,1:3) = skern;
plotkerns(1,4:5) = dkern;


int = chunkgraphinregion(cgrph,targs);

opts = []; opts.eps  = 1e-10;
sol = chunkerkerneval(cgrph,plotkerns,b,targs,opts);
u_in = skern.eval(str,targs);
u_in(int == 1) = NaN;
u_tot = u_in + sol;


figure(10);clf
h = pcolor(xx,yy,reshape(log10(abs(u_tot)),size(xx)));
h.EdgeColor = 'None';
colorbar;


%%

iverts = [1,1];
iedges = [2,1];
iedges_ref = 3-iedges;
nchs = [3,2];
cs = 0*[-1, -1];
nedge = length(cgrph.echnks);
[cgrph_fin, isort, rfins, matkerns2, plotkerns2, idignore] = chnk.fins.build_fins(cgrph, matkerns, plotkerns, iverts, iedges, iedges_ref, nchs, cs);
figure(3);clf
plot(cgrph_fin)
quiver(cgrph_fin)
hold on
for i = 1:length(rfins)
    plot(rfins{i}.r(1,:), rfins{i}.r(2,:),'.')
end
axis equal

assert(norm(cgrph_fin.r(:,:) - cgrph.r(:,isort)) == 0)

opts = []; opts.rcip_ignore = [idignore];
tic;
A2 = chunkermat(cgrph_fin,matkerns2,opts);
A2 = A2 + eye(size(A2));
toc;
opts = []; opts.rcip_ignore = [idignore]; opts.adaptive_correction = true;
tic;
A3 = chunkermat(cgrph_fin,matkerns2,opts);
A3 = A3 + eye(size(A3));
toc;


norm(A2 - A(isort,isort))
norm(A3 - A(isort,isort))
assert(norm(A3 - A(isort,isort)) < 1e-12)


B  = chunkerkerneval(cgrph,plotkerns, rhs,targs);
B2 = chunkerkerneval(cgrph_fin,plotkerns2, rhs(isort),targs);
norm(B2 - B)
assert(norm(B2 - B) < 1e-12)


% finTest1()

function finTest1()

x = linspace(-1,1,40)+1e-2; x = x(1:end-1); y = x;

[xx,yy] = meshgrid(x,y);
targs = []; targs.r = [xx(:)';yy(:)'];

verts = [0 0.3 1 1 .5;
         0 1 1 0 0];
evs = [[1;2],[2;3],[3;4],[4;5],[5;1]];
evs = flipud(evs);
cgrph = chunkgraph(verts,evs);
isplit = 2;
[cgrph_new, isort] = chnk.fins.split_chunkgraph(cgrph,isplit,3);
% figure(3);clf
% plot(cgrph_new)
% quiver(cgrph_new)


% figure(4);clf
% plot(cgrph_new.echnks(isplit))
% hold on
% plot(cgrph_new.verts(1, cgrph_new.edgesendverts(:,isplit)), cgrph_new.verts(2, cgrph_new.edgesendverts(:,isplit)), 'o')
% hold off

[cgrph_new2, isort2] = chnk.fins.split_chunkgraph(cgrph_new,1,3);


assert(norm(cgrph_new.r(:,:) - cgrph.r(:,(isort))) == 0)
assert(norm(cgrph_new2.r(:,:) - cgrph_new.r(:,isort2)) == 0)
assert(norm(cgrph_new2.r(:,:) - cgrph.r(:,isort(isort2))) == 0)

int = chunkgraphinregion(cgrph,targs);
int2 = chunkgraphinregion(cgrph_new,targs);
assert(norm(int-int2) == 0)

end