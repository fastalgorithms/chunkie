function [F] = chunkerflam(chnkr,kern,dval,opts)
%CHUNKERFLAM build the requested FLAM compressed representation 
% (e.g. a recursive skeletonization factorization) of the system matrix 
% for given kernel and chunker description of boundary. This routine
% has the same quadrature options as CHUNKERMAT.
%
% Syntax: sysmat = chunkerflam(chnkr,kern,opts)
%
% Input:
%   chnkr - chunker object describing boundary
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(srcinfo,targinfo), where 
%           srcinfo and targinfo are in the ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   dval - (default 0.0) float or float array. Let A be the matrix 
%           corresponding to on-curve convolution with the provided kernel. 
%           If a scalar is provided, the system matrix is 
%                   A + dval*eye(size(A))  
%           If a vector is provided, it should be length size(A,1). The
%           system matrix is then
%                   A + diag(dval)
%
% Optional input:
%   opts  - options structure. available options (default settings)
%           opts.flamtype = string ('rskelf'), type of compressed 
%                           representation to compute. Available:
%                               
%                               The recursive skeletonization routines
%                               should be sufficient curves which are not
%                               approximately space filling.
%
%                               - 'rskelf', recursive skeletonization
%                               factorization. Can be immediately used 
%                               for system solves (via RSKELF_SV) and 
%                               determinants (via RSKELF_LOGDET)
%                               - 'rskel', recursive skeletonization
%                               compression. Can be immediately used for
%                               matvecs or embedded in a sparse matrix.
%                               Not recommended unless you have a
%                               compelling reason.
%           opts.useproxy = boolean (true), use a proxy function to speed
%                           up the FLAM compression. It may be desirable to
%                           turn this off if there appear to be precision
%                           issues.
%
%           opts.quad = string ('ggqlog'), specify quadrature routine to 
%                       use. 
%
%                       - 'ggqlog' uses a generalized Gaussian quadrature 
%                       designed for logarithmically singular kernels and 
%                       smooth kernels with removable singularities
%                       - 'native' selects standard scaled Gauss-Legendre 
%                       quadrature for native functions
%           opts.l2scale = boolean (false), if true scale rows by 
%                           sqrt(whts) and columns by 1/sqrt(whts)
%           opts.occ = integer "occupancy" parameter (200) determines
%                   how many sources are in any leaf of the tree used to
%                   sort points. Determines the smaller dimensions of the
%                   first set of compressions. Sent to FLAM
%           opts.rank_or_tol = integer or float "rank_or_tol" 
%                   parameter (1e-14). Lower precision increases speed.
%                   Sent to FLAM
%
% Output:
%   F - the requested FLAM compressed representation of the 
%           system matrix
%
% Examples:
%   F = chunkermat(chnkr,kern,dval); % standard options
%   sysmat = chunkermat(chnkr,kern,dval,opts);
%

F = [];

if length(chnkr) > 1
    chnkr = merge(chnkr);
end

if nargin < 3
    dval = 0.0;
end

if nargin < 4
    opts = [];
end

quad = 'ggqlog';
flamtype = 'rskelf';
useproxy = true;
occ = 200;
rank_or_tol = 1e-14;


if or(chnkr.nch < 1,chnkr.k < 1)
    warning('empty chunker, doing nothing')
    sp = [];
    return
end

% determine operator dimensions using first two points

srcinfo = []; targinfo = [];
srcinfo.r = chnkr.r(:,1); srcinfo.d = chnkr.d(:,1); 
srcinfo.d2 = chnkr.d2(:,1); srcinfo.n = chnkr.n(:,1);
i2 = min(2,chnkr.npt);
targinfo.r = chnkr.r(:,i2); targinfo.d = chnkr.d(:,i2); 
targinfo.d2 = chnkr.d2(:,i2); targinfo.n = chnkr.n(:,i2);

ftemp = kern(srcinfo,targinfo);
opdims = size(ftemp);
assert(opdims(1) == opdims(2), 'the opdim should be a square matrix')

if (length(dval) == 1)
    dval = dval*ones(opdims(1)*chnkr.npt,1);
end

if (length(dval) ~= opdims(1)*chnkr.npt)
    warning('provided dval array is length %d. must be scalar or length %d',...
        length(dval),opdims(1)*chnkr.npt);
    return
end

% get opts from struct if available

if isfield(opts,'quad')
    quad = opts.quad;
end
if isfield(opts,'l2scale')
    l2scale = opts.l2scale;
end

if isfield(opts,'occ')
    occ = opts.occ;
end
if isfield(opts,'rank_or_tol')
    rank_or_tol = opts.rank_or_tol;
end

% check if chosen FLAM routine implemented before doing real work

if ~ (strcmpi(flamtype,'rskelf') || strcmpi(flamtype,'rskel'))
    warning('selected flamtype %s not available, doing nothing',flamtype)
    return
end

% get nonsmooth quadrature

if strcmpi(quad,'ggqlog') 
    chunkermatopt = struct('quad','ggq','type','log','nonsmoothonly',true);
elseif strcmpi(quad,'native')
    chunkermatopt = struct('quad','native','nonsmoothonly',true);
else
    warning('specified quadrature method not available');
    return;
end
sp = chunkermat(chnkr,kern,chunkermatopt);

[m,n] = size(sp);
sp = sp + spdiags(dval,0,m,n);

% prep and call flam

wts = weights(chnkr);
% TODO: the xflam should be repeated... 
% xflam = chnkr.r(:,:);
xflam = zeros(2,chnkr.npt*opdims(2));
for i=1:opdims(2)
    xflam(:,i:opdims(2):end) = chnkr.r(:,:);
end

width = max(max(chnkr)-min(chnkr));
optsnpxy = []; optsnpxy.rank_or_tol = rank_or_tol;
optsnpxy.nsrc = occ;

npxy = chnk.flam.nproxy_square(kern,width,optsnpxy);

matfun = @(i,j) chnk.flam.kernbyindex(i,j,chnkr,wts,kern,opdims,sp);
[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy);

if strcmpi(flamtype,'rskelf')
    ifaddtrans = true;
    pxyfun = @(x,slf,nbr,l,ctr) chnk.flam.proxyfun(slf,nbr,l,ctr,chnkr,wts, ...
        kern,opdims,pr,ptau,pw,pin,ifaddtrans);
    F = rskelf(matfun,xflam,occ,rank_or_tol,pxyfun);
end

if strcmpi(flamtype,'rskel') 
    pxyfunr = @(rc,rx,cx,slf,nbr,l,ctr) chnk.flam.proxyfunr(rc,rx,slf,nbr,l, ...
        ctr,chnkr,wts,kern,opdims,pr,ptau,pw,pin);
    F = rskel(matfun,xflam,xflam,occ,rank_or_tol,pxyfunr);
end
	 


end
