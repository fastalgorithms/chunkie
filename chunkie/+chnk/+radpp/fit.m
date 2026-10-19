function info = fit(f,rmin,rmax,opts)
%CHNK.RADPP.FIT adaptive piecewise Chebyshev fit of a vector-valued f(r).
%
% Syntax:
%   info = chnk.radpp.fit(f,rmin,rmax)
%   info = chnk.radpp.fit(f,rmin,rmax,opts)
%
% f maps a column r (nr,1) to (nr,m), possibly complex. The fit is on
% rmin < r < rmax.
%
% opts:
%   norder  - terms per panel (default 6)
%   eps     - tolerance (default 1e-12)
%   nlevmax - max levels (default ceil(log2(rmax/rmin)) + 20), warns
%   maxsub  - max panels (default 400000), errors
%   noise   - absolute error noise(j)/r allowed in component j (default 0)
%   stopnoise - stop at the noise level of f, see below (default true)
%
% A panel whose error stops halving when split, below 1e-4, is assumed
% to have stalled due to noise. Set stopnoise = false for f with sharp
% features.
%
% info: nbin, nlev, err, breaks (nbin+1,1), coefs (nbin,norder,m), the
% polynomial in r-breaks(i) on panel i, highest degree first
%
% see also CHNK.RADPP.FITEVAL, CHNK.RADPP.ROTATE

if nargin < 4 || isempty(opts), opts = struct(); end

if ~(rmin > 0) || ~(rmax > rmin)
    error('CHNK.RADPP.FIT: need 0 < rmin < rmax.');
end

norder  = 6;
tol     = 1e-12;
nlevmax = ceil(log2(rmax/rmin)) + 20;
maxsub  = 400000;
noise   = 0;
stopnoise = true;
if isfield(opts,'norder'),  norder  = opts.norder;  end
if isfield(opts,'eps'),     tol     = opts.eps;     end
if isfield(opts,'nlevmax'), nlevmax = opts.nlevmax; end
if isfield(opts,'maxsub'),  maxsub  = opts.maxsub;  end
if isfield(opts,'noise'),   noise   = opts.noise;   end
if isfield(opts,'stopnoise'), stopnoise = opts.stopnoise; end

n  = norder;
x  = cos((2*(1:n).'-1)*pi/(2*n));          % Chebyshev nodes on [-1,1]
x2 = [-1 + (x+1)/2; (x+1)/2];               % nodes of the two halves
amat = interpmat(x,x2);

% panels being refined: endpoints a, b and values fv (node, panel, comp)
a  = rmin;
b  = rmax;
fv = f(a + (x+1)*(b-a)/2);
m  = size(fv,2);
fv = reshape(fv,n,1,m);

ad = zeros(0,1);  bd = zeros(0,1);  fd = zeros(n,0,m);
err = 0;

% error of the parent of each panel
ep = inf;
inoise = false(0,1);

for ilev = 1:nlevmax

    na = numel(a);
    t  = a.' + (x2+1).*(b-a).'/2;
    fo = reshape(f(t(:)),2*n,na,m);                       % on the halves
    fi = reshape(amat*reshape(fv,n,na*m),2*n,na,m);        % interpolated

    e = relerr(fo-fi, fo, noise, tol, a);
    stag = stopnoise & e > tol & e > ep/2 & e < 1e-4;
    done = e <= tol | stag;
    if ilev == nlevmax
        done(:) = true;
        if any(e > tol)
            bad = e > tol;
            warning('CHNK:radpp:fit:nlevmax', ['CHNK.RADPP.FIT: %d panels not resolved after ' ...
                'nlevmax = %d levels, error %.1e, in r = [%.2e, %.2e].'], ...
                nnz(bad), nlevmax, max(e), min(a(bad)), max(b(bad)));
        end
    end
    err = max([err; e(done)]);

    inoise = [inoise; stag(done)];        %#ok<AGROW>
    ad = [ad; a(done)];                   %#ok<AGROW>
    bd = [bd; b(done)];                   %#ok<AGROW>
    fd = cat(2,fd,fv(:,done,:));          %#ok<AGROW>

    isplit = find(~done);
    if isempty(isplit), break; end

    if numel(ad) + 2*numel(isplit) > maxsub
        error('CHNK.RADPP.FIT: maxsub = %d panels exceeded, %d unresolved in r = [%.2e, %.2e].', ...
              maxsub, numel(isplit), min(a(isplit)), max(b(isplit)));
    end

    amid = (a(isplit) + b(isplit))/2;
    ep   = [e(isplit); e(isplit)];
    fv   = cat(2,fo(1:n,isplit,:),fo(n+1:2*n,isplit,:));
    a    = [a(isplit); amid];
    b    = [amid; b(isplit)];

end

if any(inoise)
    warning('CHNK:radpp:fit:noise', ['CHNK.RADPP.FIT: %d panels stopped at the ' ...
        'noise level of f, error %.1e > eps = %.1e, in r = [%.2e, %.2e].'], ...
        nnz(inoise), err, tol, min(ad(inoise)), max(bd(inoise)));
end

[ad,isort] = sort(ad);
bd = bd(isort);
fd = fd(:,isort,:);

nbin   = numel(ad);
breaks = [ad; bd(end)];

% monomial coefficients in r-a on each panel, highest degree first
tt = (x+1)/2;
cf = (tt.^(n-1:-1:0)) \ reshape(fd,n,nbin*m);
cf = reshape(cf,n,nbin,m) .* ((1./(bd-ad)).'.^((n-1:-1:0).'));
coefs = permute(cf,[2 1 3]);                    % (nbin,n,m)

info = struct('nbin',nbin,'nlev',ilev,'err',err, ...
              'breaks',breaks,'coefs',coefs);

end


function e = relerr(d, fo, noise, tol, a)
% largest error over nodes and components of each panel, relative to
% max(local size,1) plus the noise allowance, (npanel,1)
derr = max(max(abs(real(d)),[],1), max(abs(imag(d)),[],1));      % (1,na,m)
den  = max(max(max(abs(real(fo)),[],1), max(abs(imag(fo)),[],1)), 1);
if any(noise(:) ~= 0)
    den = den + reshape(noise(:),1,1,[])./(tol*a(:).');
end
e = reshape(max(derr./den,[],3),[],1);
end


function amat = interpmat(x,x2)
% barycentric interpolation from Chebyshev nodes x to points x2
n  = numel(x);
bw = ((-1).^((1:n).'-1)).*sin((2*(1:n).'-1)*pi/(2*n));
dd = x2 - x.';
w  = bw.'./dd;
amat = w./sum(w,2);
[i0,j0] = find(dd == 0);
for k = 1:numel(i0)
    amat(i0(k),:)     = 0;
    amat(i0(k),j0(k)) = 1;
end
end
