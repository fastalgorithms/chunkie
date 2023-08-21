
clearvars; close all;
addpaths_loc();
cparams = [];
cparams.eps = 1.0e-4;
pref = []; 
pref.k = 16;
narms = 3;
amp = 0.5;
cparams.maxchunklen = 0.1;
cparams.ta = -pi/2;
cparams.tb = pi/2;
start = tic; chnkr = chunkerfunc(@(t) starfish(t,narms,amp),cparams,pref); 
t1 = toc(start);

[~,~,info] = sortinfo(chnkr);

np = 100;
isstart = ceil(chnkr.npt/3);
isend = isstart + np;
isind = isstart:isend;


itstart = ceil(2*chnkr.npt/3);
itstart = isstart+np+1;
itend = itstart + np;
itind = itstart:itend;



srcinfo = [];
srcinfo.r = chnkr.r(:,isind);
srcinfo.d = chnkr.d(:,isind);
srcinfo.d2 = chnkr.d2(:,isind);
srcinfo.n = chnkr.n(:,isind);

targinfo = [];
targinfo.r = chnkr.r(:,itind);
targinfo.d = chnkr.d(:,itind);
targinfo.d2 = chnkr.d2(:,itind);
targinfo.n = chnkr.n(:,itind);

% Real k tests
zk = 1.1;

type = 'dprimediff';
start = tic; submat = chnk.axissymhelm2d.kern(zk, srcinfo, targinfo, type); 
t1 = toc(start);

submat2 = chnk.axissymhelm2d.kern(zk, srcinfo, targinfo, 'dprime'); 
submat2 = submat2 - chnk.axissymhelm2d.kern(1j*zk, srcinfo, targinfo, 'dprime'); 


fprintf('First call to axissymhelm2d.kern time%d\n',t1);

v = get_exact_kernels(zk, srcinfo, targinfo, type);


err1 = norm(v(:) - submat(:));
fprintf('Error in kernel %s  = %d \n',type,err1);
fprintf('ratios = %d\n', max(abs(v(1:end)./submat(1:end)-1)));
fprintf('ratios difference= %d\n', max(abs(v(1:end)./submat2(1:end)-1)));








function [v] = get_exact_kernels(zk, srcinfo, targinfo, type)

src = srcinfo.r;
targ = targinfo.r; 

srcnorm = srcinfo.n;
targnorm = targinfo.n;

[~, ns] = size(src);
[~, nt] = size(targ);

v = zeros(nt,ns);

if strcmpi(type,'s')
    fker = @(x, s, t, rns, rnt) fslp(x, zk, s, t, rns, rnt);
elseif strcmpi(type, 'd')
    fker = @(x, s, t, rns, rnt) fdlp(x, zk, s, t, rns, rnt);
elseif strcmpi(type, 'sprime')
    fker = @(x, s, t, rns, rnt) fsprime(x, zk, s, t, rns, rnt);
elseif strcmpi(type, 'dprime')
    fker = @(x, s, t, rns, rnt) fdprime(x, zk, s, t, rns, rnt);
elseif strcmpi(type, 'sdiff')
    fker = @(x, s, t, rns, rnt) fsdiff(x, zk, s, t, rns, rnt);    
elseif strcmpi(type, 'dprimediff')
    fker = @(x, s, t, rns, rnt) fdprimediff(x, zk, s, t, rns, rnt);
end

for j=1:ns
    for i=1:nt
        v(i,j) = integral(@(x) fker(x, src(:,j), targ(:,i), srcnorm(:,j), targnorm(:,i)), 0, 2*pi,'AbsTol',1E-12); 
    end
end


end


function f = fslp (x, zk, s, t, rns, rnt)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    r = sqrt(rs.^2 + rt.^2 - 2*rs.*rt.*cos(x) + (zs-zt).^2);
    f = exp(1j*zk*r)/4/pi./r.*rs;
end



function f = fdlp (x, zk, s, t, rns, rnt)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    rnd = (rt.*cos(x) - rs).*rns(1) + (zt - zs).*rns(2);
    
    r = sqrt(rs.^2 + rt.^2 - 2*rs.*rt.*cos(x) + (zs-zt).^2);
    f = rnd.*(1 - 1j*zk*r).*exp(1j*zk*r)/4/pi./r.^3.*rs;
end



function f = fsprime (x, zk, s, t, rns, rnt)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    rnd = (rt - rs.*cos(x)).*rnt(1) + (zt - zs).*rnt(2);
    
    r = sqrt(rs.^2 + rt.^2 - 2*rs.*rt.*cos(x) + (zs-zt).^2);
    f = -rnd.*(1-1j*zk*r).*exp(1j*zk*r)/4/pi./r.^3.*rs;
end



function f = fdprime (x, zk, s, t, rns, rnt)
    rs = s(1); zs = s(2);
    rt = t(1); zt = t(2);
    
    rndt = (rt - rs.*cos(x)).*rnt(1) + (zt - zs).*rnt(2);
    rnds = (rt.*cos(x) - rs).*rns(1) + (zt - zs).*rns(2);
    rnsnt = rns(1)*rnt(1).*cos(x) + rns(2)*rnt(2);
    
    r = sqrt(rs.^2 + rt.^2 - 2*rs.*rt.*cos(x) + (zs-zt).^2);
    f = -(rnsnt.*(1j*zk.*r-1).*exp(1j*zk*r)/4/pi./r.^3 + ...
           rndt.*rnds.*(-zk^2.*r.^2 - 3*1j*zk.*r + 3).*exp(1j*zk*r)/4/pi./r.^5).*rs;
end


function f = fsdiff (x, zk, s, t, rns, rnt)
    f = fslp(x, zk, s, t, rns, rnt) - fslp(x, 1j*zk, s, t, rns, rnt);
end



function f = fdprimediff (x, zk, s, t, rns, rnt)
    f1 = fdprime(x, zk, s, t, rns, rnt); 
    f2 = fdprime(x, 1j*zk, s, t, rns, rnt);
    f = f1 - f2;
end
