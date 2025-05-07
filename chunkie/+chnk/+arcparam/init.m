function param_data = init(chnkr,ich)
% CHNK.ARCPARAM.INIT initialize arclength parameterization of a chunker
% outputs param_data, which contains precomputed parameterization data
%
% Input:
%   chnkr - chunker to build parameterization
%   ich - set of panels to be parameterized. Default is the whole chunker 
%           (1:chnkr.nch). 
%
%   param_data.plen - arclength of selected panels. The arclength
%       parameterization will be defined on [0,sum(param_data.plen)]
%   param_data.pstrt - arclength coordinate of left endpoints of selected
%           panels
%   param_data.cr - arclength legendre series of chhnkr.r
%   param_data.cd - arclength legendre series of chhnkr.d
%   param_data.cd2 - arclength legendre series of chhnkr.d2
%   param_data.k - chunker order
%   param_data.dim - chunker dimension
%   param_data.nch - number of chunks in ich
%   param_data.eps - estimate of the error in the parameterization (the
%       size of largest last legendre coefficient)
%   param_data.maxcond - worst condition number encountered when
%       reparameterizing
%
% see also CHNK.ARCPARAM.EVAL

% author: Tristan Goodwill

if nargin < 2
    ich = 1:chnkr.nch;
end

% get lengths
wts = chnkr.wts(:,ich);
plen = sum(wts,1); plen = plen(:).';
pstrt = [0,cumsum(plen)];

% fetch data
rs  = chnkr.r(:,:,ich);
ds  = chnkr.d(:,:,ich);
d2s = chnkr.d2(:,:,ich);

% get arclength deriveitves
A = lege.intmat(chnkr.k);
darc = squeeze(sqrt(sum(ds.^2,1)));
s = 2*(A*darc./plen) - 1;

darc = reshape(darc,1,chnkr.k,[]);

ddarc = sum(ds.*d2s,1) ./ darc.^1;
d2s = d2s./darc.^2 - ds./darc.^3 .*ddarc;
ds = ds./darc;


% build legendre series
nchs = length(ich);

cr = zeros(chnkr.k,chnkr.dim,nchs);
cd = zeros(chnkr.k,chnkr.dim,nchs);
cd2 = zeros(chnkr.k,chnkr.dim,nchs);

maxcond=0;
for i = 1:nchs
    coef2val = lege.pols(s(:,i),chnkr.k-1).';
    [L,U] = lu(coef2val);
    cr(:,:,i)  = U\( L \ (rs(:,:,i).'));
    cd(:,:,i)  = U\( L \ (ds(:,:,i).'));
    cd2(:,:,i) = U\( L \ (d2s(:,:,i).'));
    maxcond = max(maxcond,cond(coef2val));
end

param_data = [];
param_data.plen = plen;
param_data.pstrt = pstrt;
param_data.cr = cr;
param_data.cd = cd;
param_data.cd2 = cd2;

param_data.k = chnkr.k;
param_data.dim = chnkr.dim;
param_data.nch = nchs;

% build error metrics
last_coefs = cr(end-1:end,:,:);
eps = max(last_coefs(:));
param_data.eps = eps;
param_data.maxcond = maxcond;
end

