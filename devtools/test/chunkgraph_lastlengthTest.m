% chunkgraph_lastlengthTest0();

% function chunkgraph_lastlengthTest0()

narms = 5;
rots = exp(1i*2*pi*(1:narms)/narms);
verts = zeros(2,2*narms);
radi = 1;
rado = 5;
for i = 1:narms
    verts(:,2*i-1) = radi*[real(rots(i));imag(rots(i))];
    verts(:,2*i) = rado*[real(rots(i));imag(rots(i))];
end
nverts = size(verts,2);
edge2verts = [1:nverts;circshift(1:nverts,-1)];
edge2verts = [edge2verts, [1;5]];
% curve parameters
amp = -0.5;
frq = 5;

% define curve for each edge
fchnks = cell(size(edge2verts,2)+1,1);
for i = 1:narms
    % odd edges are straight
    fchnks{2*i-1} = [];
    % even edges are curved
    fchnks{2*i} = @(t) sinearc(t,amp,frq);
end
fchnks{size(edge2verts,2)} = [];

% add extra curve
ctr = [8;0];
% verts = [verts, ctr];
edge2verts = [edge2verts, [NaN; NaN]];

cparams = cell(1,size(edge2verts,2));
for i = 1:length(cparams)
    cparams{i}.eps = 1e-8;
end

narms = 3;
amp = 0.3;
fchnks{end} = @(t) starfish(t,narms,amp, ctr, 0, 0.3);
cparams{end}.ta = 0;
cparams{end}.tb = 2*pi;

cgrph = chunkgraph(verts,edge2verts,fchnks,cparams);

plot(cgrph,'.')

ref_opts = cell(size(edge2verts,2),1);
ref_opts{1}.nover = 2;
cgrph2 = refine(cgrph,ref_opts);

plot(cgrph2,'o')

cgrph3 = cgrph;
ichnk = 2;
cgrph3.echnks(ichnk) = split(cgrph3.echnks(ichnk),3,struct('frac',2/12));

clf
plot(cgrph3,'x')



clf; hold on
last_len = 0.9;  last_len = sum(cgrph.echnks(1).wts(:,1)); last_len = 2;
cgrph4 = refine(cgrph,[],last_len);
% cgrph5 = refine(cgrph4,[],last_len);
% assert(cgrph5.npt == cgrph4.npt)
% cgrph4 = cgrph;
plot(cgrph4,'.');axis equal

hold on
% %%
cg = cgrph4;
nverts = size(verts,2);


for j = 1:nverts
% for j = 1
    loc_edges = cg.vstruc{j}{1};
    loc_dir = cg.vstruc{j}{2};

    nloc = length(loc_edges);
    if nloc == 0, continue, end
    arcs = zeros(1,nloc);

    idch_loc = [cg.echnks(loc_edges).nch];
    idch_loc(loc_dir == -1) = 1;

    for k = 1:nloc
        wts = cg.echnks(loc_edges(k)).wts;
        arcs(k) = sum(wts(:,idch_loc(k)));
        plot(cg.echnks(loc_edges(k)).r(1,:,idch_loc(k)), cg.echnks(loc_edges(k)).r(2,:,idch_loc(k)),'o')
    end
    % [diff(arcs), log2(arcs(1) / last_len), log2(arcs(2) / last_len)-round(log2(arcs(1) / last_len))]
    assert(log2(arcs(2) / last_len)-round(log2(arcs(1) / last_len)) < 1e-10)
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
