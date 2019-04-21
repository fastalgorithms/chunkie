
function chnkr = split(chnkr,ich,opts,x,w,u)
%SPLIT this routine takes the list of all chunks and splits one in
%      half with respect to arclength.
%   Input
%        ich - the chunk number to split

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


k = chnkr.k;
nch = chnkr.nch;

if nargin < 6
    [x, w, u, ~] = lege.exps(k);
end

%  first construct dsdt

dsdt = sqrt(chnkr.d(1,:,ich).^2+chnkr.d(2,:,ich).^2)*chnkr.h(ich);
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

% new points in parameter space

ts1 = -1+(t1+1)*(x+1)/2.0;
ts2 = t1+(1-t1)*(x+1)/2.0;

% evaluate the new values of r, d, d2 and 
% update nch, adj, h

xy = chnkr.r(:,:,ich);
dxy = chnkr.d(:,:,ich);
d2xy = chnkr.d2(:,:,ich);

cxy = u*(xy.');
cdxy = u*(dxy.');
cd2xy = u*(d2xy.');

i1=chnkr.adj(1,ich);
i2=chnkr.adj(2,ich);

ndim = 2;

r_1 = zeros(ndim,k,1);
d_1 = zeros(ndim,k,1);
d2_1 = zeros(ndim,k,1);
r_2 = zeros(ndim,k,1);
d_2 = zeros(ndim,k,1);
d2_2 = zeros(ndim,k,1);

for i = 1:ndim
    x = lege.exev(ts1,cxy(:,i));
    dx = lege.exev(ts1,cdxy(:,i));
    d2x = lege.exev(ts1,cd2xy(:,i));
    r_1(i,:,1) = x;
    d_1(i,:,1) = dx;
    d2_1(i,:,1) = d2x;
    x = lege.exev(ts2,cxy(:,i));
    dx = lege.exev(ts2,cdxy(:,i));
    d2x = lege.exev(ts2,cd2xy(:,i));
    r_2(i,:,1) = x;
    d_2(i,:,1) = dx;
    d2_2(i,:,1) = d2x;
end

hold=chnkr.h(ich);

% update chnkr

chnkr = chnkr.addchunk();

chnkr.h(ich) = hold*(t1+1)/2;
chnkr.h(nch+1) = hold*(1-t1)/2;

chnkr.r(:,:,ich) = r_1;
chnkr.r(:,:,nch+1) = r_2;
chnkr.d(:,:,ich) = d_1;
chnkr.d(:,:,nch+1) = d_2;
chnkr.d2(:,:,ich) = d2_1;
chnkr.d2(:,:,nch+1) = d2_2;

chnkr.adj(2,ich)=nch+1;
chnkr.adj(1,nch+1)=ich;
chnkr.adj(2,nch+1)=i2;
chnkr.adj(1,i2)=nch+1;

end
