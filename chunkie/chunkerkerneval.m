function fints = chunkerkerneval(chnkr,kern,dens,targs,opts)
%CHUNKERKERNEVAL compute the convolution of the integral kernel with
% the density defined on the chunk geometry. 
%
% Syntax: fints = chunkerkerneval(chnkr,kern,dens,targs,opts)
%
% Input:
%   chnkr - chunker object description of curve
%   kern - integral kernel taking inputs kern(srcinfo,targinfo) 
%   dens - density on boundary, should have size opdims(2) x k x nch
%          where k = chnkr.k, nch = chnkr.nch, where opdims is the 
%           size of kern for a single src,targ pair
%   targs - targ(1:2,i) gives the coords of the ith target
%
% Optional input:
%   opts - structure for setting various parameters
%       opts.flam - if = true, use flam utilities (true)
%       opts.forcesmooth - if = true, only use the smooth integration rule
%                           (false)
%       opts.forceadap - if = true, only use adaptive quadrature (false)
%    NOTE: only one of forcesmooth or forceadap is allowed. If both 
%           false, a hybrid algorithm is used, where the smooth rule is 
%           applied for targets separated by opts.fac*length of chunk from
%           a given chunk and adaptive integration is used otherwise
%       opts.fac = the factor times the chunk length used to decide 
%               between adaptive/smooth rule
%       opts.quadgkparams - if non-empty this is a cell structure
%           containing string,value pairs to be sent to quadgk (default {})
%       opts.eps = tolerance for adaptive integration
%
% output:
%   fints - opdims(1) x nt array of integral values where opdims is the 
%           size of kern for a single input (dimension of operator output)
%
% see also QUADGK

% TODO: find a method for using proxy surfaces while still not performing
%   the add/subtract strategy...

% author: Travis Askham (askhamwhat@gmail.com)

% determine operator dimensions using first two points


srcinfo = []; targinfo = [];
srcinfo.r = chnkr.r(:,1); srcinfo.d = chnkr.d(:,1); 
srcinfo.n = chnkr.n(:,1); srcinfo.d2 = chnkr.d2(:,1);
targinfo.r = chnkr.r(:,2); targinfo.d = chnkr.d(:,2); 
targinfo.d2 = chnkr.d2(:,2); targinfo.n = chnkr.n(:,2);

ftemp = kern(srcinfo,targinfo);
opdims = size(ftemp);

if nargin < 5
    opts = [];
end

forcesmooth = false;
forceadap = false;
forcepquad = false;
flam = true;
fac = 1.0;
quadgkparams = {};
eps = 1e-12;
if isfield(opts,'forcesmooth'); forcesmooth = opts.forcesmooth; end
if isfield(opts,'forceadap'); forceadap = opts.forceadap; end
if isfield(opts,'forcepquad'); forcepquad = opts.forcepquad; end
if isfield(opts,'flam'); flam = opts.flam; end
if isfield(opts,'quadgkparams'); quadgkparams = opts.quadgkparams; end
if isfield(opts,'fac'); fac = opts.fac; end
if isfield(opts,'eps'); eps = opts.eps; end

[dim,~] = size(targs);

if (dim ~= 2); warning('only dimension two tested'); end

optssmooth = []; optssmooth.flam = flam;
optsadap = []; optsadap.quadgkparams = quadgkparams;
optsadap.eps = eps;

if forcesmooth
    fints = chunkerkerneval_smooth(chnkr,kern,opdims,dens,targs, ...
        [],opts);
    return
end

if forceadap
    fints = chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
        targs,[],optsadap);
    return
end


if forcepquad
    optsflag = []; optsflag.fac = fac;
    flag = flagnear(chnkr,targs,optsflag);
    fints = chunkerkerneval_smooth(chnkr,kern,opdims,dens,targs, ...
        flag,opts);

    fints = fints + chunkerkerneval_ho(chnkr,kern,opdims,dens, ...
        targs,flag,optsadap);

    return
end

% smooth for sufficiently far, adaptive otherwise

optsflag = []; optsflag.fac = fac;
flag = flagnear(chnkr,targs,optsflag);

fints = chunkerkerneval_smooth(chnkr,kern,opdims,dens,targs, ...
    flag,opts);

fints = fints + chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
        targs,flag,optsadap);


end

function fints = chunkerkerneval_ho(chnkr,kern,opdims,dens, ...
    targs,flag,opts)

% target
[~,nt] = size(targs);
fints = zeros(opdims(1)*nt,1);
k = size(chnkr.r,2);
[t,w] = lege.exps(2*k);
ct = lege.exps(k);
bw = lege.barywts(k);
r = chnkr.r;
d = chnkr.d;
n = chnkr.n;
d2 = chnkr.d2;
h = chnkr.h;

% interpolation matrix
intp = lege.matrin(k,t);          % interpolation from k to 2*k
intp_ab = lege.matrin(k,[-1;1]);  % interpolation from k to end points
targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
targn = zeros(chnkr.dim,nt);
for j=1:size(chnkr.r,3)
    [ji] = find(flag(:,j));
    if(~isempty(ji))
        idxjmat = (j-1)*k+(1:k);


        % Helsing-Ojala (interior/exterior?)
        mat1 = chnk.pquadwts(r,d,n,d2,h,ct,bw,j,targs(:,ji), ...
            targd(:,ji),targn(:,ji),targd2(:,ji),kern,opdims,t,w,opts,intp_ab,intp); % depends on kern, different mat1?

        fints(ji) = fints(ji) + mat1*dens(idxjmat);
    end
end
end


function fints = chunkerkerneval_smooth(chnkr,kern,opdims,dens, ...
    targs,flag,opts)

flam = true;

if nargin < 6
    flag = [];
end
if nargin < 7
    opts = [];
end
if isfield(opts,'flam'); flam = opts.flam; end

k = chnkr.k;
nch = chnkr.nch;

assert(numel(dens) == opdims(2)*k*nch,'dens not of appropriate size')
dens = reshape(dens,opdims(2),k,nch);

[~,w] = lege.exps(k);
[~,nt] = size(targs);

fints = zeros(opdims(1)*nt,1);

targinfo = []; targinfo.r = targs;

% assume smooth weights are good enough

if ~flam
    % do dense version
    if isempty(flag)
        % nothing to ignore
        for i = 1:nch
            densvals = dens(:,:,i); densvals = densvals(:);
            dsdtdt = sqrt(sum(abs(chnkr.d(:,:,i)).^2,1));
            dsdtdt = dsdtdt(:).*w(:)*chnkr.h(i);
            dsdtdt = repmat( (dsdtdt(:)).',opdims(2),1);
            densvals = densvals.*(dsdtdt(:));
            srcinfo = []; srcinfo.r = chnkr.r(:,:,i); 
            srcinfo.n = chnkr.n(:,:,i);
            srcinfo.d = chnkr.d(:,:,i); srcinfo.d2 = chnkr.d2(:,:,i);
            kernmat = kern(srcinfo,targinfo);
            fints = fints + kernmat*densvals;
        end
    else
        % ignore interactions in flag array
        for i = 1:nch
            densvals = dens(:,:,i); densvals = densvals(:);
            dsdtdt = sqrt(sum(abs(chnkr.d(:,:,i)).^2,1));
            dsdtdt = dsdtdt(:).*w(:)*chnkr.h(i);
            dsdtdt = repmat( (dsdtdt(:)).',opdims(2),1);
            densvals = densvals.*(dsdtdt(:));
            srcinfo = []; srcinfo.r = chnkr.r(:,:,i); 
            srcinfo.n = chnkr.n(:,:,i);
            srcinfo.d = chnkr.d(:,:,i); srcinfo.d2 = chnkr.d2(:,:,i);
            kernmat = kern(srcinfo,targinfo);

            rowkill = find(flag(:,i)); 
            rowkill = (opdims(1)*(rowkill(:)-1)).' + (1:opdims(1)).';
            kernmat(rowkill,:) = 0;

            fints = fints + kernmat*densvals;
        end
    end
else

%     % hack to ignore interactions: create sparse array with very small
%     % number instead of zero in ignored entries. kernbyindexr overwrites
%     %   with that small number
%     [i,j] = find(flag);
%     
%     % targets correspond to multiple outputs (if matrix valued kernel)
%     inew = (opdims(1)*(i(:)-1)).' + (1:opdims(1)).'; inew = inew(:);
%     jnew = (j(:)).' + 0*(1:opdims(1)).'; jnew = jnew(:);
%     
%     % source chunks have multiple points and can be multi-dimensional
%     jnew = (opdims(2)*chnkr.k*(jnew(:)-1)).' + (1:(chnkr.k*opdims(2))).';
%     jnew = jnew(:);
%     inew = (inew(:)).' + 0*(1:(chnkr.k*opdims(2))).';
%     inew = inew(:);
%     
%     mm = nt*opdims(1); nn = chnkr.npt*opdims(2);
%     v = 1e-300*ones(length(inew),1); sp = sparse(inew,jnew,v,mm,nn);

    wts = weights(chnkr);
    wts = wts(:);
    
    xflam1 = chnkr.r(:,:);
    xflam1 = repmat(xflam1,opdims(2),1);
    xflam1 = reshape(xflam1,chnkr.dim,numel(xflam1)/chnkr.dim);
    targsflam = repmat(targs(:,:),opdims(1),1);
    targsflam = reshape(targsflam,chnkr.dim,numel(targsflam)/chnkr.dim);
%    matfun = @(i,j) chnk.flam.kernbyindexr(i,j,targs,chnkr,wts,kern, ...
%        opdims,sp);
    matfun = @(i,j) chnk.flam.kernbyindexr(i,j,targs,chnkr,wts,kern, ...
        opdims);
    

    width = max(abs(max(chnkr)-min(chnkr)))/3;
    tmax = max(targs(:,:),[],2); tmin = min(targs(:,:),[],2);
    wmax = max(abs(tmax-tmin));
    width = max(width,wmax/3);
    npxy = chnk.flam.nproxy_square(kern,width);
    [pr,ptau,pw,pin] = chnk.flam.proxy_square_pts(npxy);

    pxyfun = @(rc,rx,cx,slf,nbr,l,ctr) chnk.flam.proxyfunr(rc,rx,slf,nbr,l, ...
        ctr,chnkr,wts,kern,opdims,pr,ptau,pw,pin);
    

    optsifmm=[]; optsifmm.Tmax=Inf;
    F = ifmm(matfun,targsflam,xflam1,200,1e-14,pxyfun,optsifmm);
    fints = ifmm_mv(F,dens(:),matfun);

    % delete interactions in flag array (possibly unstable approach)
    
    targinfo = [];
    if ~isempty(flag)
        for i = 1:nch
            densvals = dens(:,:,i); densvals = densvals(:);
            dsdtdt = sqrt(sum(abs(chnkr.d(:,:,i)).^2,1));
            dsdtdt = dsdtdt(:).*w(:)*chnkr.h(i);
            dsdtdt = repmat( (dsdtdt(:)).',opdims(2),1);
            densvals = densvals.*(dsdtdt(:));
            srcinfo = []; srcinfo.r = chnkr.r(:,:,i); 
            srcinfo.n = chnkr.n(:,:,i);
            srcinfo.d = chnkr.d(:,:,i); srcinfo.d2 = chnkr.d2(:,:,i);

            delsmooth = find(flag(:,i)); 
            delsmoothrow = (opdims(1)*(delsmooth(:)-1)).' + (1:opdims(1)).';
            delsmoothrow = delsmoothrow(:);
            targinfo.r = targs(:,delsmooth);

            kernmat = kern(srcinfo,targinfo);
            fints(delsmoothrow) = fints(delsmoothrow) - kernmat*densvals;
        end
    end    
end

end

function fints = chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
    targs,flag,opts)

k = chnkr.k;
nch = chnkr.nch;

if nargin < 6
    flag = [];
end
if nargin < 7
    opts = [];
end

assert(numel(dens) == opdims(2)*k*nch,'dens not of appropriate size')
dens = reshape(dens,opdims(2),k,nch);

[~,nt] = size(targs);

fints = zeros(opdims(1)*nt,1);

% using adaptive quadrature

[t,w] = lege.exps(2*k+1);
ct = lege.exps(k);
bw = lege.barywts(k);
r = chnkr.r;
d = chnkr.d;
n = chnkr.n;
d2 = chnkr.d2;
h = chnkr.h;
targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
targn = zeros(chnkr.dim,nt);

if isempty(flag) % do all to all adaptive
    for i = 1:nch
        for j = 1:nt
            fints1 = chnk.adapgausskerneval(r,d,n,d2,h,ct,bw,i,dens,targs(:,j), ...
                    targd(:,j),targn(:,j),targd2(:,j),kern,opdims,t,w,opts);
            
            indj = (j-1)*opdims(1);
            ind = indj(:).' + (1:opdims(1)).'; ind = ind(:);
            fints(ind) = fints(ind) + fints1;
        end
    end
else % do only those flagged
    for i = 1:nch
        [ji] = find(flag(:,i));
        [fints1,maxrec,numint,iers] =  chnk.adapgausskerneval(r,d,n,d2,h,ct,bw,i,dens,targs(:,ji), ...
                    targd(:,ji),targn(:,ji),targd2(:,ji),kern,opdims,t,w,opts);
                
        indji = (ji-1)*opdims(1);
        ind = (indji(:)).' + (1:opdims(1)).';
        ind = ind(:);
        fints(ind) = fints(ind) + fints1;        

    end
    
end

end
