% chunkgraph_lastlengthTest0();

% function chunkgraph_lastlengthTest0()

narms = 3;
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

% curve parameters
amp = -0.5;
frq = 5;

% define curve for each edge
fchnks = cell(2*narms,1);
for i = 1:narms
    % odd edges are straight
    fchnks{2*i-1} = [];
    % even edges are curved
    fchnks{2*i} = @(t) sinearc(t,amp,frq);
end

cgrph = chunkgraph(verts,edge2verts,fchnks);

plot(cgrph,'.')

ref_opts = cell(2*narms,1);
ref_opts{1}.nover = 2;
cgrph2 = refine(cgrph,ref_opts);

plot(cgrph2,'o')

cgrph3 = cgrph;
ichnk = 2;
cgrph3.echnks(ichnk) = split(cgrph3.echnks(ichnk),3,struct('frac',2/12));

clf
plot(cgrph3,'x')



clf
last_len = 0.9;  last_len = sum(cgrph.echnks(1).wts(:,1));
cgrph4 = refine(cgrph,[],last_len);
cgrph5 = refine(cgrph4,[],last_len);
assert(cgrph5.npt == cgrph4.npt)
% cgrph4 = cgrph;
plot(cgrph4,'.')
hold on
% %%
% cg = cgrph4;
% nverts = narms*2;
% for j = 1:nverts
%     loc_edges = cg.vstruc{j}{1};
%     loc_dir = cg.vstruc{j}{2};
% 
%     nloc = length(loc_edges);
%     arcs = zeros(1,nloc);
% 
%     idch_loc = [cg.echnks(loc_edges).nch];
%     idch_loc(loc_dir == 1) = 1;
% 
%     for k = 1:nloc
%         wts = cg.echnks(loc_edges(k)).wts;
%         arcs(k) = sum(wts(:,idch_loc(k)));
%         plot(cg.echnks(loc_edges(k)).r(1,:,idch_loc(k)), cg.echnks(loc_edges(k)).r(2,:,idch_loc(k)),'o')
%     end
%     [diff(arcs), log2(arcs(1) / last_len), log2(arcs(1) / last_len)-round(log2(arcs(1) / last_len))]
% end



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
