chunkgraph_basicTest0();
dyadic_chunkgraphTest0()


function chunkgraph_basicTest0()
%chunkgraph_basicTest
%
% test of basic chunkgraph class constructors, methods, etc.
%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   .  .  .  builds a simple pentagonal chunkergraph 


verts = exp(1i*2*pi*(0:4)/5);
verts = [real(verts);imag(verts)];

endverts = [1:5; [2:5 1]];

edge2verts = [-1, 1, 0, 0, 0; ...
               0,-1, 1, 0, 0; ...
               0, 0,-1, 1, 0; ...
               0, 0, 0,-1, 1; ...
               1, 0, 0, 0,-1];

edge2verts = sparse(edge2verts);
amp = 0.5;
frq = 6;
fchnks    = {};
for icurve = 1:size(edge2verts,1)
    fchnks{icurve} = @(t) sinearc(t,amp,frq);
end

[cgrph1] = chunkgraph(verts,edge2verts,fchnks);
[cgrph2] = chunkgraph(verts,endverts,fchnks);

assert(nnz(cgrph1.v2emat-cgrph2.v2emat) == 0);

cgrph1 = balance(cgrph1);

% a graph with two edges in the old format

verts = [1 0 1; -1 0 1]; edge2verts = [-1 1 0; 0 -1 1];
cg1 = chunkgraph(verts,edge2verts);

edgesendverts = [1:2; 2:3];
cg2 = chunkgraph(verts,edgesendverts);

assert(nnz(cg1.v2emat-cg2.v2emat) == 0);

% testing multiply connected

verts = [1 0 -1 2 0 -2; 0 1 0 -0.5 2 -0.5]; 
edgends = [1 2 3 4 5 6; 2 3 1 5 6 4];
cg1 = chunkgraph(verts,edgends);

assert(numel(cg1.regions) == 3);

% testing connected by bridge

verts = [1 0 -1 4 3 2; 0 1 0 0 1 0]; 
edgends = [1 2 3 4 5 6 1; 3 1 2 6 4 5 6];
cg1 = chunkgraph(verts,edgends);

assert(numel(cg1.regions) == 3)

% testing a loop

verts = [2;1]; edgends = [1;1];
cfun = @(t) [cos(t(:).'); sin(t(:).').*sin(0.5*t(:).')];
fchnks = {cfun};
cparams = []; cparams.ta = 0; cparams.tb = 2*pi;
pref = []; pref.k = 12;
cg1 = chunkgraph(verts,edgends,fchnks,cparams);

assert(numel(cg1.regions) == 2)

% testing nested 

verts = [1 0 -1 2 0 -2; 0 1 0 -1 2 -1];
edgends = [1 2 3 4 5 6; 2 3 1 5 6 4];
cg1 = chunkgraph(verts,edgends);

assert(numel(cg1.regions) == 3)

% region id for targets: compare to inpolygon routine for polygonal domains

% adjacent triangles

verts = [1 0 -1 0; 0 1 0 -1]; edgesendverts = [1:3, 3, 4; 2:3, 1, 4, 1];
cg = chunkgraph(verts,edgesendverts);

x1 = linspace(-pi,pi); 
[xx,yy] = meshgrid(x1,x1);

targs = [xx(:).'; yy(:).'];

ids = chunkgraphinregion(cg,targs);
ids2 = chunkgraphinregion(cg,{x1,x1});

idstrue = polygonids(cg,xx,yy);

assert(nnz(ids(:)-idstrue(:)) == 0)
assert(nnz(ids(:)-ids2(:)) == 0)

A = [3 2; 1 1]; 
v = [-1; 2];
cg = A*cg + v;
targs = A*targs + v;

ids = chunkgraphinregion(cg,targs);

assert(nnz(ids(:)-idstrue(:)) == 0)

cg = cg*2;
targs = 2*targs;

ids = chunkgraphinregion(cg,targs);

assert(nnz(ids(:)-idstrue(:)) == 0)

cg = cg.rotate(pi/4);
targs = [cos(pi/4) -sin(pi/4); sin(pi/4) cos(pi/4)]*targs;

ids = chunkgraphinregion(cg,targs);
assert(nnz(ids(:)-idstrue(:)) == 0)

cg = cg.reflect(pi/2);
targs(1,:) = -targs(1,:);

ids = chunkgraphinregion(cg,targs);
assert(nnz(ids(:)-idstrue(:)) == 0)

% nested triangles

verts = [1 0 -1 2 0 -2; 0 1 0 -1 2 -1];
edgends = [1 2 3 4 5 6; 2 3 1 5 6 4];
cg = chunkgraph(verts,edgends);

x1 = linspace(-pi,pi); % avoid edge cases...
[xx,yy] = meshgrid(x1,x1);

targs = [xx(:).'; yy(:).'];

ids = chunkgraphinregion(cg,targs);

idstrue = polygonids(cg,xx,yy);

assert(nnz(ids(:)-idstrue(:)) == 0)

%% build a chunkgraph and refine 

% adjacent triangles

verts = [1 0 -1 0; 0 1 0 -1]; edgesendverts = [1:3, 3, 4; 2:3, 1, 4, 1];
cg = chunkgraph(verts,edgesendverts);
refopts = []; refopts.nover = 1;
nchs1 = [cg.echnks.nch];
cg = refine(cg,refopts);
nchs2 = [cg.echnks.nch];

assert(all(2*nchs1 == nchs2))

%% 

verts = exp(1i*2*pi*(0:4)/5);
verts = [real(verts);imag(verts)];

endverts = [1:5; [2:5 1]];

amp = 0.5;
frq = 6;
fchnks = @(t) sinearc(t,amp,frq);
% simpler calling sequence if all edges same function
cgrph = chunkgraph(verts,endverts,fchnks);



end

function dyadic_chunkgraphTest0()
% test option to dyadically refine chunkgraph into specific corners

n = 5;
thetas = (1:n)/n*2*pi;
verts = [cos(thetas); sin(thetas)];
edge2vert = [1:n; circshift(1:n,-1)];

% vertices to refine
vert_ref = [1,4];

% init cparams
cparams = cell(1,n);
for i = 1:n
    cparams{i}.chsmall = [Inf,Inf];
end

edgeids = 1:n;
% find edges that start at a vertex that will be refined
edge_ref1 = edgeids(ismember(edge2vert(1,:),vert_ref));
for i = edge_ref1
    cparams{i}.chsmall(1) = 1e-8;
end

% find edges that end at a vertex that will be refined
edge_ref2 = edgeids(ismember(edge2vert(2,:),vert_ref));
for i = edge_ref2
    cparams{i}.chsmall(2) = 1e-8;
end

cgrph = chunkgraph(verts,edge2vert,[],cparams);

lens = chunklen(cgrph.echnks(edge_ref1(1)));
assert(lens(1) < cparams{edge_ref1(1)}.chsmall(1))
lens = chunklen(cgrph.echnks(edge_ref2(2)));
assert(lens(end) < cparams{edge_ref2(2)}.chsmall(2))

end


function idstrue = polygonids(cg,xx,yy)

verts = cg.verts;
edgesendverts = cg.edgesendverts;

idstrue = nan(size(xx(:)));
regions = cg.regions;
for j = 1:numel(regions)
    vlist = [];
    for ic = 1:numel(regions{j})
        elist = regions{j}{ic};
        for ie = 1:numel(elist)
            eid = elist(ie);
            if eid > 0
                v1 = edgesendverts(1,eid);
                v2 = edgesendverts(2,eid);
            else
                v1 = edgesendverts(2,-eid);
                v2 = edgesendverts(1,-eid);
            end
            vlist = [vlist, v1];
            if ie == numel(elist)
                vlist = [vlist, v2];
            end
        end
    end
    intmp = inpolygon(xx(:),yy(:),verts(1,vlist),verts(2,vlist));
    if j == 1
        idstrue(intmp == 0) = j;
    else
        idstrue(intmp > 0) = j;
    end
end

end

function [r,d,d2] = sinearc(t,amp,frq)

xs = t;
ys = amp*sin(frq*t);
xp = ones(size(t));
yp = amp*frq*cos(frq*t);
xpp = zeros(size(t));
ypp = -frq*frq*amp*sin(t);

r = [(xs(:)).'; (ys(:)).'];
d = [(xp(:)).'; (yp(:)).'];
d2 = [(xpp(:)).'; (ypp(:)).'];
end
