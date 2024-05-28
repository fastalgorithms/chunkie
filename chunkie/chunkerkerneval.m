function fints = chunkerkerneval(chnkobj,kern,dens,targobj,opts)
%CHUNKERKERNEVAL compute the convolution of the integral kernel with
% the density defined on the chunk geometry. 
%
% Syntax: fints = chunkerkerneval(chnkr,kern,dens,targs,opts)
%
% Input:
%   chnkobj - chunker object or chunkgraph object description of curve
%   kern - kernel class object or kernel function taking inputs 
%                      kern(srcinfo,targinfo) 
%   dens - density on boundary, should have size opdims(2) x k x nch
%          where k = chnkr.k, nch = chnkr.nch, where opdims is the 
%           size of kern for a single src,targ pair
%   targobj - object describing the target points, can be specified as
%       * array of points
%       * chunker object
%       * chunkgraph object
%
% Optional input:
%   opts - structure for setting various parameters
%       opts.flam - if = true, use flam utilities. to be replaced by the 
%                   opts.forceflam flag. 
%                   opts.flam supercedes opts.accel, if
%                   both are true, then flam will be used. (false)
%       opts.accel - if = true, use specialized fmm if defined 
%                   for the kernel, if it doesnt exist or if too few 
%                   sources/targets, or if false, 
%                   do direct. (true)
%       opts.forcefmm - if = true, use specialized fmm if defined,
%                   independent of the number of sources/targets. (false)
%       opts.forcesmooth - if = true, only use the smooth integration rule
%                           (false)
%       opts.forceadap - if = true, only use adaptive quadrature (false)
%    NOTE: only one of forcesmooth or forceadap is allowed. If both 
%           false, a hybrid algorithm is used, where the smooth rule is 
%           applied for targets separated by opts.fac*length of chunk from
%           a given chunk and adaptive integration is used otherwise
%       opts.fac = the factor times the chunk length used to decide 
%               between adaptive/smooth rule
%       opts.eps = tolerance for adaptive integration
%
% output:
%   fints - opdims(1) x nt array of integral values where opdims is the 
%           size of kern for a single input (dimension of operator output)
%

% TODO: find a method for using proxy surfaces while still not performing
%   the add/subtract strategy...

% author: Travis Askham (askhamwhat@gmail.com)

% determine operator dimensions using first two points

if isa(kern,'kernel')
    kerneval = kern.eval;
else
    kerneval = kern;
end

% Assign appropriate object to chnkr
if class(chnkobj) == "chunker"
   chnkr = chnkobj;
elseif class(chnkobj) == "chunkgraph"
   chnkr = merge(chnkobj.echnks);
else
    msg = "Unsupported object in chunkerkerneval";
    error(msg)
end

srcinfo = []; targinfo = [];
srcinfo.r = chnkr.r(:,1); srcinfo.d = chnkr.d(:,1); 
srcinfo.n = chnkr.n(:,1); srcinfo.d2 = chnkr.d2(:,1);
targinfo.r = chnkr.r(:,2); targinfo.d = chnkr.d(:,2); 
targinfo.d2 = chnkr.d2(:,2); targinfo.n = chnkr.n(:,2);

ftemp = kerneval(srcinfo,targinfo);
opdims = size(ftemp);

if nargin < 5
    opts = [];
end

opts_use = [];
opts_use.forcesmooth = false;
opts_use.forceadap = false;
opts_use.forcepquad = false;
opts_use.flam = false;
opts_use.accel = true;
opts_use.forcefmm = false;
opts_use.fac = 1.0;
opts_use.eps = 1e-12;
if isfield(opts,'forcesmooth'); opts_use.forcesmooth = opts.forcesmooth; end
if isfield(opts,'forceadap'); opts_use.forceadap = opts.forceadap; end
if isfield(opts,'forcepquad'); opts_use.forcepquad = opts.forcepquad; end
if isfield(opts,'flam')
    opts_use.flam = opts.flam;
end
if isfield(opts,'accel'); opts_use.accel = opts.accel; end
if isfield(opts,'fac'); opts_use.fac = opts.fac; end
if isfield(opts,'eps'); opts_use.eps = opts.eps; end

% Assign appropriate object to targinfo
targinfo = [];
if isa(targobj, "chunker")
    targinfo.r = targobj.r(:,:);
    targinfo.d = targobj.d(:,:);
    targinfo.d2 = targobj.d2(:,:);
    targinfo.n = targobj.n(:,:);
elseif isa(targobj, "chunkgraph")
    targinfo.r = targobj.r(:,:);
    targinfo.d = targobj.d(:,:);
    targinfo.d2 = targobj.d2(:,:);
    targinfo.n = targobj.n(:,:);
else
    targinfo.r = targobj;
end



[dim,~] = size(targinfo.r);

if (dim ~= 2); warning('only dimension two tested'); end

if opts_use.forcesmooth
    fints = chnk.chunkerkerneval_smooth(chnkr,kern,opdims,dens,targinfo, ...
        [],opts_use);
    return
end

if opts_use.forceadap
    fints = chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
        targinfo,[],opts_use);
    return
end


if opts_use.forcepquad
    optsflag = []; optsflag.fac = opts_use.fac;
    flag = flagnear(chnkr,targinfo.r,optsflag);
    fints = chnk.chunkerkerneval_smooth(chnkr,kern,opdims,dens,targinfo, ...
        flag,opts_use);

    fints = fints + chunkerkerneval_ho(chnkr,kern,opdims,dens, ...
        targinfo,flag,opts_use);

    return
end

% smooth for sufficiently far, adaptive otherwise

rho = 1.8;
optsflag = [];  optsflag.rho = rho;
flag = flagnear_rectangle(chnkr,targinfo.r,optsflag);

npoly = chnkr.k*2;
nlegnew = chnk.ellipse_oversample(rho,npoly,opts_use.eps);
nlegnew = max(nlegnew,chnkr.k);

[chnkr2,dens2] = upsample(chnkr,nlegnew,dens);
fints = chnk.chunkerkerneval_smooth(chnkr2,kern,opdims,dens2,targinfo, ...
    flag,opts_use);

fints = fints + chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
        targinfo,flag,opts_use);

end

function fints = chunkerkerneval_ho(chnkr,kern,opdims,dens, ...
    targinfo,flag,opts)

% target
[~,nt] = size(targinfo.r);
fints = zeros(opdims(1)*nt,1);
k = size(chnkr.r,2);
[t,w] = lege.exps(2*k);
ct = lege.exps(k);
bw = lege.barywts(k);
r = chnkr.r;
d = chnkr.d;
n = chnkr.n;
d2 = chnkr.d2;

% interpolation matrix
intp = lege.matrin(k,t);          % interpolation from k to 2*k
intp_ab = lege.matrin(k,[-1;1]);  % interpolation from k to end points

targs = targinfo.r;

targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
targn = zeros(chnkr.dim,nt);
if isfield(targinfo, 'd')
    targd = targinfo.d;
end

if isfield(targinfo, 'd2')
    targd2 = targinfo.d2;
end

if isfield(targinfo, 'n')
    targn = targinfo.n;
end
for j=1:size(chnkr.r,3)
    [ji] = find(flag(:,j));
    if(~isempty(ji))
        idxjmat = (j-1)*k+(1:k);


        % Helsing-Ojala (interior/exterior?)
        mat1 = chnk.pquadwts(r,d,n,d2,ct,bw,j,targs(:,ji), ...
            targd(:,ji),targn(:,ji),targd2(:,ji),kern,opdims,t,w,opts,intp_ab,intp); % depends on kern, different mat1?

        fints(ji) = fints(ji) + mat1*dens(idxjmat);
    end
end
end



function fints = chunkerkerneval_adap(chnkr,kern,opdims,dens, ...
    targinfo,flag,opts)

if isa(kern,'kernel')
    kerneval = kern.eval;
else
    kerneval = kern;
end

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

% Extract target info
targs = targinfo.r;
[~,nt] = size(targs);
targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
targn = zeros(chnkr.dim,nt);
if isfield(targinfo, 'd')
    targd = targinfo.d;
end

if isfield(targinfo, 'd2')
    targd2 = targinfo.d2;
end

if isfield(targinfo, 'n')
    targn = targinfo.n;
end


fints = zeros(opdims(1)*nt,1);

% using adaptive quadrature

[t,w] = lege.exps(2*k+1);
ct = lege.exps(k);
bw = lege.barywts(k);
r = chnkr.r;
d = chnkr.d;
n = chnkr.n;
d2 = chnkr.d2;


if isempty(flag) % do all to all adaptive
    for i = 1:nch
        for j = 1:nt
            fints1 = chnk.adapgausskerneval(r,d,n,d2,ct,bw,i,dens,targs(:,j), ...
                    targd(:,j),targn(:,j),targd2(:,j),kerneval,opdims,t,w,opts);
            
            indj = (j-1)*opdims(1);
            ind = indj(:).' + (1:opdims(1)).'; ind = ind(:);
            fints(ind) = fints(ind) + fints1;
        end
    end
else % do only those flagged
    for i = 1:nch
        [ji] = find(flag(:,i));
        [fints1,maxrec,numint,iers] =  chnk.adapgausskerneval(r,d,n,d2,ct,bw,i,dens,targs(:,ji), ...
                    targd(:,ji),targn(:,ji),targd2(:,ji),kerneval,opdims,t,w,opts);
                
        indji = (ji-1)*opdims(1);
        ind = (indji(:)).' + (1:opdims(1)).';
        ind = ind(:);
        fints(ind) = fints(ind) + fints1;        

    end
    
end

end
