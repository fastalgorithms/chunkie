function ells = ellipses(chnkr,rho,opts)

nch = chnkr.nch;
nth = max(2*nch,20);

k = chnkr.k;
nch = chnkr.nch;
dim = chnkr.dim;

assert(dim == 2,'ellipse images only well defined in 2d');

zell = lege.bernstein_ellipse(nth,rho);

rcoef = exps(chnkr);
zcoef = reshape(rcoef(1,:,:)+1i*rcoef(2,:,:),k,nch);

km1=k-1;
zpols = lege.pols(zell,km1);
zpols = zpols.';

ells = zpols*zcoef;
ells = reshape([real(ells(:).');imag(ells(:).')],2,nth,nch);
