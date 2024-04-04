%% variable opdims test

clearvars; close all;

nverts = 4; 
verts = exp(1i*2*pi*(0:(nverts-1))/nverts);
verts = [real(verts);imag(verts)];

iind = 1:nverts;
jind = 1:nverts;

iind = [iind iind];
jind = [jind jind + 1];
jind(jind>nverts) = 1;
svals = [-ones(1,nverts) ones(1,nverts)];
edge2verts = sparse(iind,jind,svals,nverts,nverts);


fchnks    = cell(1,size(edge2verts,1));

cparams = [];
cparams.nover = 2;
[cgrph] = chunkgraph(verts,edge2verts,fchnks,cparams);

vstruc = procverts(cgrph);
rgns = findregions(cgrph);
cgrph = balance(cgrph);


nregions = 2;
zk = 30;

kerns(4,4) = kernel();
kerns(1:2,1:2) =  2*kernel('helm', 'd', zk); % Dirichlet on 1,2
kerns(3:4,1:2) = -2*kernel('helm', 'dp', zk); % Neumann on 1,2


kerns(1:2,3:4) =  2*kernel('helm', 's', zk);
kerns(3:4,3:4) = -2*kernel('helm', 'sp', zk);


start = tic; sysmat = chunkermat(cgrph,kerns); t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

sys = eye(size(sysmat,1)) + sysmat;

fkernsrc = kernel('helm','s',zk);
sources = [1;1];
strengths = [1];
srcinfo = []; srcinfo.r = sources;
targinfo = []; targinfo.r = merge(cgrph.echnks(1:2)).r(:,:); 
targinfo.d = merge(cgrph.echnks(1:2)).d(:,:);
bdrydatad = fkernsrc.fmm(1e-12,srcinfo,targinfo,strengths);

fkernsrc = kernel('helm','sp',zk);
targinfo = []; targinfo.r = merge(cgrph.echnks(3:4)).r(:,:); 
targinfo.d = merge(cgrph.echnks(3:4)).d(:,:);
targinfo.n = merge(cgrph.echnks(3:4)).n(:,:);
bdrydatan = fkernsrc.fmm(1e-12,srcinfo,targinfo,strengths);
bdrydata = [bdrydatad;bdrydatan];

sol = sys\bdrydata;


% 
% udense = sys*bdry_data;
% cormat = chunkermat(cgrph,kerns,struct("corrections",true));
% sysapply = @(sigma) bieapply(cgrph,kerns,1,sigma,cormat);
% start = tic; u = sysapply(bdry_data); t1 = toc(start);
% fprintf('%5.2e s : time for matrix free apply\n',t1)
% 
% relerr = norm(udense-u)/norm(udense);
% fprintf('relative apply error %5.2e\n',relerr);
% assert(relerr < 1e-13)
% 
% sol1 = sys\bdry_data;
% sol2 = gmres(sysapply, bdry_data, [], 1e-12, 200);
% 
% relerr = norm(sol1-sol2)/norm(sol1);
% fprintf('relative solve error %5.2e\n',relerr);
% assert(relerr < 1e-10)



% nverts = 4; 
% verts = exp(1i*2*pi*(0:(nverts-1))/nverts);
% verts = [real(verts);imag(verts)];
% 
% iind = 1:nverts;
% jind = 1:nverts;
% 
% iind = [iind iind];
% jind = [jind jind + 1];
% jind(jind>nverts) = 1;
% svals = [-ones(1,nverts) ones(1,nverts)];
% edge2verts = sparse(iind,jind,svals,nverts,nverts);
% edge2verts = [edge2verts; [-1,0,1,0]];
% 
% fchnks    = cell(1,size(edge2verts,1));
% 
% cparams = [];
% cparams.nover = 2;
% [cgrph] = chunkgraph(verts,edge2verts,fchnks,cparams);
% 
% vstruc = procverts(cgrph);
% rgns = findregions(cgrph);
% cgrph = balance(cgrph);
% 
% 
% nregions = 2;
% ks = [1.1;2.1]*30;
% cs = [];
% coefs = [1;1];
% cc1 = -ones(2,2);
% cc2 = ones(2,2);
% 
% kerns(5,5) = kernel();
% kerns(1:2,1:2) = kernel('helm', 'd', ks(1)); % Dirichlet on 1,2
% kerns(3:4,1:2) = kernel('helm', 'dp', ks(2)); % Neumann on 1,2
% 
% 
% kerns(1:2,3:4) = kernel('helm', 's', ks(1));
% kerns(3:4,3:4) = kernel('helm', 'sp', ks(2));
% 
% kerntmp = @(s,t) -chnk.helm2d.kern(ks(1),s,t,'all',cc1)- ...
%                  chnk.helm2d.kern(ks(2),s,t,'all',cc2);
% kerns(5,5) = kernel(kerntmp);
