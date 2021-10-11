function submat = diagbuildmat(r,d,n,d2,h,i,fkern,opdims,...
				      glwts,logquad)
%%CHNK.QUADJH.DIAGBUILDMAT   
% kernel-split quadrature for self-interaction of logarithmic singularity
% type.
rs = r(:,:,i); ds = d(:,:,i); d2s = d2(:,:,i); hs = h(i); 
ns = n(:,:,i);

srcinfo = []; srcinfo.r = rs; srcinfo.d = ds; srcinfo.d2 = d2s; srcinfo.n = ns;
targinfo = []; targinfo.r = rs; targinfo.d = ds; targinfo.d2 = d2s; targinfo.n = ns;

dsnrms = sqrt(sum(ds.^2,1));
ws = kron(hs(:),glwts(:));
dsdt = dsnrms(:).*ws;

dsdtndim2 = repmat(dsdt(:).',opdims(2),1);
dsdtndim2 = dsdtndim2(:);

LogC = logquad.LogC;

mat = fkern(srcinfo,targinfo);

mat = bsxfun(@times,mat,(dsdtndim2(:)).');

% the following line works for Helmhotz kernels when the argument of the
% Hankel functions is real.
submat = mat-2*imag(mat).*LogC/pi;
% corrections for diagonal entries. This part is different for different
% kernels such as D, S, D', and S'.
% The following corrections assume that fkern = [D, S; D', S'];

N = length(glwts)*opdims(2);
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
submat(inds) = 0.5*(tmp-(LogC(inds)+log(hs*dsnrms))/pi).*dsdt.';

% 3. diagonal corrections for D'. Note well: this is only correct for the
% difference kernel D'_{k1}-D'_{k2}!!!
inddp = 2:2*(N+1):N^2;
tmp = 1i/2-(log(omega/2)+gamma-0.5)/pi;
submat(inddp) = omega^2*0.25*(tmp-(LogC(inddp)+log(hs*dsnrms))/pi).*dsdt.';

% 4. diagonal corrections for S'. There is a minus sign in front, then
% there is another minus sign in targnorm and srcnorm. Hence the same sign
% for the diagonal corrections for D and S'.
indsp = N+2:2*(N+1):N^2;
submat(indsp) = submat(indd);

