function submat = buildmat(chnkr,fkern,opdims,glwts,ilist,logquad)
%%CHNK.QUADJH.BUILDMAT   
% kernel-split quadrature for self-interaction of logarithmic singularity
% type.

ds = chnkr.d;
d2s = chnkr.d2;

srcinfo = []; srcinfo.r = chnkr.r; srcinfo.d = chnkr.d;
srcinfo.d2 = chnkr.d2; srcinfo.n = chnkr.n;
targinfo = srcinfo;
hs = chnkr.h;

dsnrms = sqrt(sum(chnkr.d.^2,1)); dsnrms = dsnrms(:);
ws = kron(hs(:),glwts(:));

k = chnkr.k;
nch = chnkr.nch;
hds = zeros(1,k*nch);
for i=1:nch
  hds((i-1)*k+(1:k)) = hs(i)*dsnrms((i-1)*k+(1:k));
end

dsdt = dsnrms(:).*ws;

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

if nch == 1
  LogC = logquad.LogC;
else
  if round(hs(1)/hs(2)) == 1
    LogC = logquad.LogC1;
  elseif round(hs(1)/hs(2)) == 2
    LogC = logquad.LogC0;
  end
end

mat = fkern(srcinfo,targinfo);

% compute only needed matrix blocks. somehow this is slower.

% blist = ilist(:).';
% glist=setdiff(1:nch,blist);
% mat = zeros(k*nch*opdims(2));
% srcg = []; srcg.r=chnkr.r(:,:,glist);srcg.d=chnkr.d(:,:,glist);
% srcg.d2 = chnkr.d2(:,:,glist); srcg.n = chnkr.n(:,:,glist);
% 
% srcb = []; srcb.r=chnkr.r(:,:,blist);srcb.d=chnkr.d(:,:,blist);
% srcb.d2 = chnkr.d2(:,:,blist); srcb.n = chnkr.n(:,:,blist);
% 
% ind0 = (1:k*opdims(2));
% n0 = k*opdims(2);
% indg = [];
% for i=glist
%   indg = [indg (i-1)*n0+ind0];
% end
% indb = [];
% for i=blist
%   indb = [indb (i-1)*n0+ind0];
% end
% 
% mat(indg,indg) = fkern(srcg,srcg);
% mat(indb,indg) = fkern(srcg,srcb);
% mat(indg,indb) = fkern(srcb,srcg);

submat = bsxfun(@times,mat,(dsdtndim2(:)).');

% the following line works for Helmhotz kernels when the argument of the
% Hankel functions is real.
% myind=find(LogC);
% submat(myind) = submat(myind)-2*imag(submat(myind)).*LogC(myind)/pi;
submat = submat -2*imag(submat).*LogC/pi;

% corrections for diagonal entries. This part is different for different
% kernels such as D, S, D', and S'.
% The following corrections assume that fkern = [D, S; D', S'];

N = opdims(2)*k*nch;
% 1. diagonal corrections for D
dnrm3 = dsnrms.^3;
kappa = bsxfun(@rdivide,ds(1,:).*d2s(2,:)-ds(2,:).*d2s(1,:),dnrm3(:).');
indd = 1:2*(N+1):N^2;
submat(indd) = -kappa.*dsdt.'/4/pi; 

% 2. diagonal corrections for S. Note that hs is actually half of the 
% interval length. 
omega = logquad.omega;
gamma = 0.5772156649015328606;
tmp = 1i/2-(log(omega/2)+gamma)/pi;
inds = N+1:2*(N+1):N^2;
%submat(inds) = 0.5*(tmp-(LogC(inds)+log(hs*dsnrms))/pi).*dsdt.';
submat(inds) = 0.5*(tmp-(LogC(inds)+log(hds))/pi).*dsdt.';

% 3. diagonal corrections for D'. Note well: this is only correct for the
% difference kernel D'_{k1}-D'_{k2}!!!
inddp = 2:2*(N+1):N^2;
tmp = 1i/2-(log(omega/2)+gamma-0.5)/pi;
%submat(inddp) = omega^2*0.25*(tmp-(LogC(inddp)+log(hs*dsnrms))/pi).*dsdt.';
submat(inddp) = omega^2*0.25*(tmp-(LogC(inddp)+log(hds))/pi).*dsdt.';

% 4. diagonal corrections for S'. There is a minus sign in front, then
% there is another minus sign in targnorm and srcnorm. Hence the same sign
% in the diagonal corrections for D and S'.
indsp = N+2:2*(N+1):N^2;
submat(indsp) = submat(indd);

