function [res_flag] = datares(obj,opts)
%DATARES check if the data in the specified data rows are resolved
% By default, the data in a row on a chunk are resolved if the
% 1-norm of the last chnkr.k/2 coefficients of the Legendre coefficients 
% for that data row is less than the tolerance (1e-6). 
%
%   Let c(i,j) denote the jth Legendre coefficient of the idata(i) data
%   row on chunk l. Then, setting v(i) = \|c(i,chnkr.k/2+1:end)\|_1
%   we have 
%
%   res_flag(i,l) = v(i) < 1e-6
%
% Variants of 
% most of these options are available. 
%
% Syntax: res_flag = datares(chnkr,opts)
%
% Input:
%   chnkr - chunker object
%   opts - options data structure. 
%       opts.idata - the data rows to check resolution on 
%                   (default all rows, i.e. 1:chnkr.datadim)
%       opts.ncoeff = number of final coeffs to take the norm of
%       opts.pleg = integer or Inf indicating what norm to apply to 
%                   the Legendre coefficients (default 1)
%       opts.tol = tolerance (default 1e-6)
%       opts.pscale = integer, scale tolerance by 
%                       1/(length of chunk)^pscale for each chunk
%                       (default 0)
%       opts.rel = boolean, scale the tolerance by same norm applied to
%                   the other coefficients, i.e. the first k-ncoeff (false)
% Output:
%   res_flag(i,l) - for each chunk l and data row idata(i) is true if 
%       resolved, false otherwise
%
% Examples:
%   opts = []; opts.idata = 2:3; % checks resolution of 2nd and 3rd data
%                                % rows
%   res_flag = datares(chnkr,opts);
% 
% see also REFINE

% author: Travis Askham (askhamwhat@gmail.com)

wts = weights(obj);
lens = sum(wts,1);

idata = 1:obj.datadim;
ncoeff = floor((obj.k+0.1)/2);
pleg = 1;
tol = 1e-6;
pscale = 0;
rel = false;

if nargin < 2
    opts = [];
end

if isfield(opts,'idata')
    idata = opts.idata;
end
if isfield(opts,'ncoeff')
    ncoeff = opts.ncoeff;
end
if isfield(opts,'pleg')
    pleg = opts.pleg;
end
if isfield(opts,'tol')
    tol = opts.tol;
end
if isfield(opts,'pscale')
    pscale = opts.pscale;
end
if isfield(opts,'rel')
    rel = opts.rel;
end

nch = obj.nch;
k = obj.k;

[~,~,u] = lege.exps(k);

data = obj.data(idata,:,:);
icoeff1 = (k-ncoeff+1):k;
icoeff2 = 1:(k-ncoeff);

res_flag = false(length(idata),nch);

if ~rel
    u1 = u(icoeff1,:);

    for i=1:nch
        datai = data(:,:,i);
        cf1 = u1*(datai.');
        cf1n = vecnorm(cf1,pleg,1); cf1n = cf1n(:);
        
        res_flag(:,i) = cf1n*lens(i)^pscale < tol;
    end
else
    u1 = u(icoeff1,:);
    u2 = u(icoeff2,:);

    for i=1:nch
        datai = data(:,:,i);
        cf1 = u1*(datai.');
        cf2 = u2*(datai.');
        cf1n = vecnorm(cf1,pleg,1); cf1n = cf1n(:);
        cf2n = vecnorm(cf2,pleg,1); cf2n = cf2n(:);
        
        res_flag(:,i) = cf1n*lens(i)^pscale < tol*cf2n;
    end
    
end
        