function in = chunkerin(chnkr,pts,opts)
%CHUNKERIN returns an array indicating whether each point in
% pts is inside the domain
%
% TODO: faster algorithm (smooth rule + compare normal at nearest point)

assert(chnkr.dim == 2,'interior only well-defined for 2D');

if nargin < 3
    opts = [];
end

if ~isfield(opts,'gausseps')
    opts.gausseps = 1e-8;
end
if ~isfield(opts,'quadgkparams')
    opts.quadgkparams = {};
end
if ~isfield(opts,'verb')
    opts.verb = false;
end
if ~isfield(opts,'justsmoothworks')
    opts.justsmoothworks = false;
end

kernd = @(s,t,sn,tn) glapkern(s,t,sn,tn,'d');
dens1 = ones(chnkr.k,chnkr.nch);

ndims = [1 1];

opts.usesmooth=true;
d1 = chunkerintkern(chnkr,kernd,ndims,dens1,pts,opts); 

eps = opts.gausseps;
smoothworks = or(abs(d1) < eps,abs(d1+1) < eps);

if opts.justsmoothworks
    in = smoothworks;
    return
end

ptsfail = pts(:,~smoothworks);
if opts.verb; fprintf('npts adaptive %d\n',nnz(~smoothworks)); end
opts.usesmooth=false;
d12 = chunkerintkern(chnkr,kernd,ndims,dens1,ptsfail,opts); 

d1(~smoothworks) = d12;

in = abs(d1+1) < eps;

