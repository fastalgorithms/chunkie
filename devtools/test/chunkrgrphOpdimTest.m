clear all
addpaths_loc();

verts = [-3 1; 3 1; 3 5; -3 5; 
         -3 -5; 3 -5; 3 -1; -3 -1].';

edge2verts = [
    -1, 1, 0, 0, 0, 0, 0, 0;
    0, -1, 1, 0, 0, 0, 0, 0;
    0, 0, -1, 1, 0, 0, 0, 0;
    1, 0, 0, -1, 0, 0, 0, 0;
    0, 0, 0, 0, -1, 1, 0, 0;
    0, 0, 0, 0, 0, -1, 1, 0;
    0, 0, 0, 0, 0, 0, -1, 1;
    0, 0, 0, 0, 1, 0, 0, -1];

edge2verts = sparse(edge2verts);

fchnks = {}; % this by default gives me straight lines

prefs = struct('maxchunklen',0.5);
[cgrph] = chunkgraph(verts, edge2verts, fchnks, prefs);

% figure(1); clf; hold on; axis equal; axis off;
% plot(cgrph);
% quiver(cgrph);
% hold off;

vstruc = procverts(cgrph);
rgns = findregions(cgrph);
cgrph = balance(cgrph);

% kerns

zk0 = 1; % exterior
zk1 = 3; % interior
coef = [1 -1j*real(zk0)]; 
cc = [1 1; 1 1]; % these numbers might be wrong... 

fkern11 = @(s,t) chnk.helm2d.kern(zk0,s,t,'all',cc) - chnk.helm2d.kern(zk1,s,t,'all',cc);
fkern12 = @(s,t) chnk.helm2d.kern(zk0,s,t,'c2trans',coef);
fkern21 = @(s,t) chnk.helm2d.kern(zk0,s,t,'trans_rep');
fkern22  = @(s,t) chnk.helm2d.kern(zk0,s,t,'c',[1,1i]);

chnk_trans_flag = [1 1 1 1 0 0 0 0];

fkerns = {};

ncurve = 8;

for it=1:ncurve
    for is=1:ncurve
        kern_type = [chnk_trans_flag(it) chnk_trans_flag(is)];
        if all(kern_type==[1 1]), fkerns{it,is} = fkern11; end
        if all(kern_type==[1 0]), fkerns{it,is} = fkern12; end
        if all(kern_type==[0 1]), fkerns{it,is} = fkern21; end
        if all(kern_type==[0 0]), fkerns{it,is} = fkern22; end
    end
end

opts = struct('nonsmoothonly',false, 'rcip',true);
[sysmat] = chunkermat(cgrph, fkerns, opts);