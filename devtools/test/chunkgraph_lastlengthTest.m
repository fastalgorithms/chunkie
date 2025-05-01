chunkgraph_lastlengthTest0();

function chunkgraph_lastlengthTest0()
% test the refine option to ensure that all arclengths touching each vertex 
% agree and are a negative power of 2 times last_len

% setup vertics
ncircedge = 6;
rots = exp(1i*2*pi*(1:ncircedge)/ncircedge);
verts = zeros(2,ncircedge);
rad = 6;
for i = 1:ncircedge
    verts(:,i) = rad*[real(rots(i));imag(rots(i))];
end
nverts = size(verts,2);

% setup edges
edge2verts = [1:nverts;circshift(1:nverts,-1)];
edge2verts = [edge2verts, [2;6]];

% curve parameters
amp = -0.5;
frq = 6;

% define curve for each edge
fchnks = cell(size(edge2verts,2)+1,1);
for i = 1:ncircedge/2
    % odd edges are straight
    fchnks{2*i-1} = [];
    % even edges are curved
    fchnks{2*i} = @(t) sinearc(t,amp,frq);
end
fchnks{size(edge2verts,2)} = [];

% add extra closed curve
ctr = [0;0];
% verts = [verts, ctr];
edge2verts = [edge2verts, [NaN; NaN]];
narm = 3; amp = 0.3;
fchnks{end} = @(t) starfish(t,narm,amp, ctr, 0, 0.3);

% set cparams
cparams = cell(1,size(edge2verts,2));
for i = 1:length(cparams)
    cparams{i}.eps = 1e-8;
end

cparams{end}.ta = 0;
cparams{end}.tb = 2*pi;

% build cgrph
cgrph = chunkgraph(verts,edge2verts,fchnks,cparams);

% check that we can refine a single edge
ref_opts = [];
ref_opts.nover = 1;
ref_opts.dlist = 2;
cgrph2 = refine(cgrph,ref_opts);

ref_opts = [];
ref_opts.splitchunks = cell(1,length(cgrph.echnks));
ref_opts.splitchunks{3} = 3;
cgrph3 = refine(cgrph,ref_opts);

% test refine with last_len
last_len = 2;
cgrph4 = refine(cgrph,struct("last_len",last_len));

ifplot = 0;
if ifplot
    figure(1); clf;
    plot(cgrph4,'.'); hold on
end

% verify that all arclengths touching each vertex agree and are a negative
% power of 2 times last_len
for j = 1:nverts
    loc_edges = cgrph4.vstruc{j}{1};
    loc_dir = cgrph4.vstruc{j}{2};

    nloc = length(loc_edges);
    if nloc == 0, continue, end
    arcs = zeros(1,nloc);

    idch_loc = [cgrph4.echnks(loc_edges).nch];
    idch_loc(loc_dir == -1) = 1;

    for k = 1:nloc
        wts = cgrph4.echnks(loc_edges(k)).wts;
        arcs(k) = sum(wts(:,idch_loc(k)));
        if ifplot
        plot(cgrph4.echnks(loc_edges(k)).r(1,:,idch_loc(k)), cgrph4.echnks(loc_edges(k)).r(2,:,idch_loc(k)),'o')
        end
    end
    assert(norm(diff(arcs)) < 1e-8)
    assert(log2(arcs(1) / last_len)-round(log2(arcs(1) / last_len)) < 1e-10)
    assert(round(log2(arcs(1) / last_len)) <= 0)
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
