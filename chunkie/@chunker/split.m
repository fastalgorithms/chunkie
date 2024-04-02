function chnkr = split(chnkr,ich,opts,x,w,u,stype)
%SPLIT this routine takes the list of all chunks and splits one in
%      half with respect to either arclength or parameter space
%       (this routine is not necessarily designed to be user-callable,
%       though it is fairly simple)
%
% Syntax: chnkr = split(chnkr,ich,opts,x,w,u,stype)
%
% Input:
%   chnkr - the chunker object
%   ich - the chunk number to split
%
% Optional input:
%   opts - options structure
%       opts.thresh = threshold for newton (1e-10)
%       opts.nitermax = maximum iters for newton (1000)
%   x - precomputed Legendre nodes of order chnkr.k
%   w - precomputed Legendre weights
%   u - precomputed vals at legendre nodes -> coeffs matrix
%   stype - type of split ('a'), 'a' arclength split, 
%               't' parameter space split.
%
% Output:
%   chnkr - modified chunker object
%
% Examples:
%   [x,w,u] = lege.exps(chnkr.k); % precompute relevant quantities
%   chnkr = split(chnkr,10,[],x,w,u,'a')
%   chnkr = split(chnkr,10); % routine computes x,w,u
%
% see also REFINE

% author: Travis Askham (askhamwhat@gmail.com)

stype1 = 'a';
if nargin >= 7
    stype1 = stype;
end

r = chnkr.rstor(:,:,ich);
d = chnkr.dstor(:,:,ich);
d2 = chnkr.d2stor(:,:,ich);

thresh=1.0e-10;
nitermax = 1000;

if nargin < 3
    opts = [];
end

if isfield(opts,'thresh')
    thresh = opts.thresh;
end
if isfield(opts,'nitermax')
    nitermax = opts.nitermax;
end

nch = chnkr.nch;

if nargin < 6
    [x, w, u, ~] = lege.exps(chnkr.k);
end

t1 = 0;

if strcmpi(stype1,'a')
%  first construct dsdt

dsdt = sqrt(sum(d.^2,1));
dsdt = dsdt(:);
rltot = dot(dsdt,w);

cdsdt = u*dsdt;

t1=0;
rlhalf=rltot/2;

% use Newton to find t such that length(-1,t) = length(t,1) = rl/2

ifdone = 0;
for ijk = 1:nitermax
    ts = -1+(t1+1)*(x + 1)/2.0;
    ws = (t1+1)*w/2.0;

    vals = lege.exev(ts,cdsdt);
    rl1 = dot(vals,ws);
    val = lege.exev(t1,cdsdt);
    err=rl1-rlhalf;
    if (abs(err) < thresh); ifdone=ifdone+1; end
    if (ifdone >= 3) 
        break;
    end
    t1=t1-(rl1-rlhalf)/val;
end

if (ifdone < 3); warning('did not converge'); end

end



% new points in parameter space

ts1 = -1+(t1+1)*(x+1)/2.0;
ts2 = t1+(1-t1)*(x+1)/2.0;

% evaluate the new values of r, d, d2 and 
% update nch, adj, h

cr = u*(r.');
cd = u*(d.');
cd2 = u*(d2.');

r_1 = lege.exev(ts1,cr);
r_1 = r_1.';
r_2 = lege.exev(ts2,cr);
r_2 = r_2.';
d_1 = lege.exev(ts1,cd);
d_1 = d_1.';
d_2 = lege.exev(ts2,cd);
d_2 = d_2.';
d2_1 = lege.exev(ts1,cd2);
d2_1 = d2_1.';
d2_2 = lege.exev(ts2,cd2);
d2_2 = d2_2.';

% update chnkr

%i1=chnkr.adjstor(1,ich);
i2=chnkr.adjstor(2,ich);

chnkr = chnkr.addchunk();

h1 = (t1+1)/2;
h2 = (1-t1)/2;

chnkr.rstor(:,:,ich) = r_1;
chnkr.rstor(:,:,nch+1) = r_2;
chnkr.dstor(:,:,ich) = d_1*h1;
chnkr.wtsstor(:,ich) = (sqrt(sum(d_1.^2,1)).') .* chnkr.wstor;
chnkr.dstor(:,:,nch+1) = d_2*h2;
chnkr.wtsstor(:,nch+1) = (sqrt(sum(d_2.^2,1)).') .* chnkr.wstor;
chnkr.d2stor(:,:,ich) = d2_1*h1*h1;
chnkr.d2stor(:,:,nch+1) = d2_2*h2*h2;

chnkr.adjstor(2,ich)=nch+1;
chnkr.adjstor(1,nch+1)=ich;
chnkr.adjstor(2,nch+1)=i2;
if i2 > 0
    chnkr.adjstor(1,i2)=nch+1;
end

if (i2 < 0 && numel(chnkr.vert)~=0)
    ii = (chnkr.vert{-i2} == ich);
    chnkr.vert{-i2}(ii) = nch+1;
end

if chnkr.hasdata
    data = chnkr.datastor(:,:,ich);
    cdata = u*(data.');
    data1 = lege.exev(ts1,cdata);
    data1 = data1.';
    data2 = lege.exev(ts2,cdata);
    data2 = data2.';
    chnkr.datastor(:,:,ich) = data1;
    chnkr.datastor(:,:,nch+1) = data2;
end


end
