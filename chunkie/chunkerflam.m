function [F] = chunkerflam(chnkobj,kern,dval,opts)
%CHUNKERFLAM build the requested FLAM compressed representation 
% (e.g. a recursive skeletonization factorization) of the system matrix 
% for given kernel and chunker description of boundary. This routine
% has the same quadrature options as CHUNKERMAT.
%
% Syntax: sysmat = chunkerflam(chnkobj,kern,dval,opts)
%
% Input:
%   chnkobj - chunker object or chunkgraph object describing boundary 
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
%           opts.verb = boolean (false), if true print out the process of
%                   the compression of FLAM.
%           opts.lvlmax = integer (inf), maximum level of compression
%           opts.adaptive_correction, boolean (false), see chunkermat
% 
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

if isa(kern,'function_handle')
    kern2 = kernel(kern);
    kern = kern2;
elseif isa(kern,'cell')
    sz = size(kern);
    kern2(sz(1),sz(2)) = kernel();
    for j = 1:sz(2)
        for i = 1:sz(1)
            if isa(kern{i,j},'function_handle')
                kern2(i,j) = kernel(kern{i,j});
            elseif isa(kern{i,j},'kernel')
                kern2(i,j) = kern{i,j};
            else
                msg = "Second input is not a kernel object, function handle, " ...
                    + "or cell array";
                error(msg);
            end
        end
    end
    kern = kern2;
    
elseif ~isa(kern,'kernel')
    msg = "Second input is not a kernel object, function handle, " ...
                + "or cell array";
    error(msg);
end

if nargin < 3
    dval = 0.0;
end

if nargin < 4
    opts = [];
end

quad = 'ggqlog';        if isfield(opts,'quad'),    quad = opts.quad;end
l2scale = false;        if isfield(opts,'l2scale'), l2scale = opts.l2scale;end
flamtype = 'rskelf';    if isfield(opts,'flamtype'),flamtype = opts.flamtype;end
useproxy = true;        if isfield(opts,'useproxy'),useproxy = opts.useproxy;end
occ = 200;              if isfield(opts,'occ'),     occ = opts.occ;end
rank_or_tol = 1e-14;    if isfield(opts,'rank_or_tol'),rank_or_tol = opts.rank_or_tol;end
verb = false;           if isfield(opts,'verb'),    verb = opts.verb;end
lvlmax = inf;           if isfield(opts,'lvlmax'),  lvlmax = opts.lvlmax;end
adaptive_correction = false; if isfield(opts,'adaptive_correction'), ...
    adaptive_correction = opts.adaptive_correction;end

eps = rem(rank_or_tol,1);


% Flag for determining whether input object is a chunkergraph
icgrph = 0;

if (class(chnkobj) == "chunker")
    chnkrs = chnkobj;
elseif(class(chnkobj) == "chunkgraph")
    icgrph = 1;
    chnkrs = chnkobj.echnks;
else
    msg = "Unsupported object in chunkermat";
    error(msg)
end

for chnkr = chnkrs
    if or(chnkr.nch < 1,chnkr.k < 1)
        warning('empty chunker, doing nothing')
        return
    end
end

nchunkers = length(chnkrs);
opdims_mat = zeros(2,nchunkers,nchunkers);
lchunks    = zeros(nchunkers,1);

for i=1:nchunkers
    
    targinfo = [];
   	targinfo.r = chnkrs(i).r(:,2); targinfo.d = chnkrs(i).d(:,2); 
   	targinfo.d2 = chnkrs(i).d2(:,2); targinfo.n = chnkrs(i).n(:,2);
    lchunks(i) = size(chnkrs(i).r(:,:),2);
    
    for j=1:nchunkers
        
        % determine operator dimensions using first two points

        srcinfo = []; 
        srcinfo.r = chnkrs(j).r(:,1); srcinfo.d = chnkrs(j).d(:,1); 
        srcinfo.d2 = chnkrs(j).d2(:,1); srcinfo.n = chnkrs(j).n(:,1);

        
        if (size(kern) == 1)
            ftemp = kern.eval(srcinfo,targinfo);
        else
            ktmp = kern(i,j).eval;
            ftemp = ktmp(srcinfo,targinfo);
        end   
        opdims = size(ftemp);
        opdims_mat(:,i,j) = opdims;
    end
end    

irowlocs = zeros(nchunkers+1,1);
icollocs = zeros(nchunkers+1,1);

irowlocs(1) = 1;
icollocs(1) = 1;
for i=1:nchunkers
   icollocs(i+1) = icollocs(i) + lchunks(i)*opdims_mat(2,1,i);
   irowlocs(i+1) = irowlocs(i) + lchunks(i)*opdims_mat(1,i,1);
end    

nrows = irowlocs(end)-1;
ncols = icollocs(end)-1;

assert(nrows==ncols,'chunkerflam: nrows and ncols must be equal');

if (length(dval) == 1)
    dval = dval*ones(nrows,1);
elseif (length(dval) ~= nrows)
    warning('provided dval array is length %d. must be scalar or length %d',...
        length(dval),nrows);
    return
end

% check if chosen FLAM routine implemented before doing real work

if ~ (strcmpi(flamtype,'rskelf') || strcmpi(flamtype,'rskel'))
    warning('selected flamtype %s not available, doing nothing',flamtype)
    return
end

% get nonsmooth quadrature

if strcmpi(quad,'ggqlog') 
    chunkermatopt = struct('quad','ggq','type','log','nonsmoothonly',true, ...
        'l2scale',l2scale,'adaptive_correction',adaptive_correction,'eps',eps);
elseif strcmpi(quad,'native')
    chunkermatopt = struct('quad','native','nonsmoothonly',true, ...
        'l2scale',l2scale,'adaptive_correction',adaptive_correction,'eps',eps);
else
    warning('specified quadrature method not available');
    return;
end
sp = chunkermat(chnkobj,kern,chunkermatopt);
sp = sp + spdiags(dval,0,nrows,nrows);

% prep and call flam

wts = zeros(ncols,1);
xflam = zeros(2,ncols);

mmax = [-inf;-inf];
mmin = [inf;inf];

for i=1:nchunkers
    chnkr = chnkrs(i);
    opdim = opdims_mat(2,1,i);
    xi = chnkr.r(:,:);
    for j=1:opdim
        xflam(:,icollocs(i)+j-1:opdim:icollocs(i+1)-1) = xi;
    end
    mmax = max([mmax,max(chnkr)],[],2);
    mmin = min([mmin,min(chnkr)],[],2);
end

width = max(mmax-mmin);

chnkrtotal = merge(chnkrs);

matfun = @(i,j) chnk.flam.kernbyindex(i,j,chnkrs,kern,opdims_mat,sp,l2scale);

if ~useproxy
    if strcmpi(flamtype,'rskelf')
        F = rskelf(matfun,xflam,occ,rank_or_tol,[],struct('verb',verb,'lvlmax',lvlmax));
    end
    if strcmpi(flamtype,'rskel') 
        F = rskel(matfun,xflam,xflam,occ,rank_or_tol,[]);
    end    
    return
end

optsnpxy = []; optsnpxy.rank_or_tol = rank_or_tol;
optsnpxy.nsrc = occ;

npxy = chnk.flam.nproxy_square(kern,width,optsnpxy);

if npxy == -1
    warning('chunkerflam: proxy failed, defaulting to no proxy')
    if strcmpi(flamtype,'rskelf')
        F = rskelf(matfun,xflam,occ,rank_or_tol,[],struct('verb',verb,'lvlmax',lvlmax));
    end
    if strcmpi(flamtype,'rskel') 
        F = rskel(matfun,xflam,xflam,occ,rank_or_tol,[]);
    end    
    return
end


[pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy);

if strcmpi(flamtype,'rskelf')
    ifaddtrans = true;
    pxyfun = @(x,slf,nbr,l,ctr) chnk.flam.proxyfun(slf,nbr,l,ctr,chnkrs, ...
        kern,opdims_mat,pr,ptau,pw,pin,ifaddtrans,l2scale);
    F = rskelf(matfun,xflam,occ,rank_or_tol,pxyfun,struct('verb',verb,...
        'lvlmax',lvlmax));
end

if strcmpi(flamtype,'rskel') 
    warning('chunkerflam: proxyr has not adapted for the multiple chunker case yet')
    pxyfunr = @(rc,rx,cx,slf,nbr,l,ctr) chnk.flam.proxyfunr(rc,rx,slf,nbr,l, ...
        ctr,chnkr,kern,opdims,pr,ptau,pw,pin);
    F = rskel(matfun,xflam,xflam,occ,rank_or_tol,pxyfunr);
end
	 


end
