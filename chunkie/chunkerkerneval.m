function fints = chunkerkerneval(chnkr,kern,dens,targs,opts)
%CHUNKERKERNEVAL compute the convolution of the integral kernel with
% the density defined on the chunk geometry. 
%
% input:
%   chnkr - chunks description of curve
%   kern - integral kernel taking inputs kern(srcinfo,targinfo) 
%   opdims - input and output dimensions of the kernel (opdims(1) dimension
%           output, opdims(2) dimension of input)
%   dens - density on boundary, should have size opdims(2) x k x nch
%          where k = chnkr.k, nch = chnkr.nch
%   targs - targ(1:2,i) gives the coords of the ith target
%   opts - structure for setting various parameters
%       opts.usesmooth - if = 1, then just use the smooth integration
%          rule for each chunk. if = 0, adaptive integration 
%          (quadgk) is used. if = 2, a hybrid method is used 
%          where the smooth rule is used for points where the 
%          smooth rule is accurate to opts.gausseps digits for 
%          gauss's id and uses adaptive at other points (default 0)
%       opts.quadgkparams - if non-empty this is a cell structure
%       containing string,value pairs to be sent to quadgk (default {})
%       opts.gausseps - if the hybrid method is used, the smooth 
%       rule is applied at points for which Gauss' ID is accurate 
%       with the smooth rule to absolute error opts.gausseps (default 1e-8)
%
% output:
%   fints - opdims(1) x nt array of integral values
%
% see also QUADGK

% determine operator dimensions using first two points


srcinfo = []; targinfo = [];
srcinfo.r = chnkr.r(:,1); srcinfo.d = chnkr.d(:,1); 
srcinfo.d2 = chnkr.d2(:,1);
targinfo.r = chnkr.r(:,2); targinfo.d = chnkr.d(:,2); 
targinfo.d2 = chnkr.d2(:,2);

ftemp = kern(srcinfo,targinfo);
opdims = size(ftemp);

if nargin < 6
    opts = [];
end

[dim,nt] = size(targs);
assert(dim==2,'only dimension two tested');

if ~isfield(opts,'usesmooth'); opts.usesmooth = false; end
if ~isfield(opts,'quadgkparams'); opts.quadgkparams = {}; end
if ~isfield(opts,'gausseps'); opts.gausseps = 1e-8; end
if ~isfield(opts,'verb'); opts.verb= false; end

if opts.usesmooth == 1
    fints = chunkerkerneval_smooth(chnkr,kern,opdims,dens, ...
        targs);
elseif opts.usesmooth == 0
    fints = chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
        targs,opts);
elseif opts.usesmooth == 2
    fints = zeros(opdims(1),nt);
    optssw = []; optssw.gausseps = opts.gausseps; 
    optssw.justsmoothworks = true;
    sw = chunkerinterior(chnkr,targs,optssw);
    fints(:,sw) = reshape(...
        chunkerkerneval_smooth(chnkr,kern,opdims,dens, ...
        targs(:,sw)),opdims(1),nnz(sw));
    fints(:,~sw) = reshape(...
        chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
        targs(:,~sw),opts),opdims(1),nnz(~sw));
    fints = fints(:);
end

end



function fints = chunkerkerneval_smooth(chnkr,kern,opdims,dens, ...
    targs)

k = chnkr.k;
nch = chnkr.nch;

assert(numel(dens) == opdims(2)*k*nch,'dens not of appropriate size')
dens = reshape(dens,opdims(2),k,nch);

[~,w] = lege.exps(k);
[~,nt] = size(targs);

fints = zeros(opdims(1)*nt,1);

targinfo = []; targinfo.r = targs;

% assume smooth weights are good enough
for i = 1:nch
    densvals = dens(:,:,i); densvals = densvals(:);
    dsdtdt = sqrt(sum(abs(chnkr.d(:,:,i)).^2,1));
    dsdtdt = dsdtdt(:).*w(:)*chnkr.h(i);
    dsdtdt = repmat( (dsdtdt(:)).',opdims(2),1);
    densvals = densvals.*(dsdtdt(:));
    srcinfo = []; srcinfo.r = chnkr.r(:,:,i); 
    srcinfo.d = chnkr.d(:,:,i); srcinfo.d2 = chnkr.d2(:,:,i);
    kernmat = kern(srcinfo,targinfo);

    fints = fints + kernmat*densvals;
end

end

function fints = chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
    targs,opts)

k = chnkr.k;
nch = chnkr.nch;

assert(numel(dens) == opdims(2)*k*nch,'dens not of appropriate size')
dens = reshape(dens,opdims(2),k,nch);

[~,~,u] = lege.exps(k);
[~,nt] = size(targs);

fints = zeros(opdims(1)*nt,1);

% using adaptive quadrature
[rc,dc,d2c] = exps(chnkr);
for i = 1:nch
    if opts.verb; fprintf('chunk %d integral\n',i); end
    rci = rc(:,:,i);
    dci = dc(:,:,i);
    d2ci = d2c(:,:,i);    
    densvals = dens(:,:,i); densvals = densvals.';
    densc = u*densvals; % each column is set of coefficients
                    % for one dimension of density on chunk
    for j = 1:nt
        indj = (j-1)*opdims(1);
        for l = 1:opdims(1)
            ind = indj+l;
            temp = chunkerintchunk_kernfcoefs(kern,opdims,l,...
                densc,rci,dci,d2ci,targs(:,j));

            fints(ind) = fints(ind) + temp*chnkr.h(i);
        end
    end
end

end
