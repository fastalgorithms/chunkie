function [chnkr,eps] = arc_param_chunker(chnkr,opts)
% ARC_PARAM_CHUNKER reparameterize chnkr by arclength
%
% NOTE: requires an accurate chunker discretization
% Input:
%   chnkr - chunker to be reparameterized
%   opts - struct of opts (default option)
%       opts.mv_bdries = boolean (false). If true, panel boundaries are 
%           allowed to move. maxchunklen is garaunteed not to increase
%
% Output: reparameterized chnkr and tolerance (eps) used to construct it
%
% see also CHUNKERARCPARAM

% author: Tristan Goodwill

if nargin == 1
    opts = [];
end

if ~isfield(opts,'mv_bdries')
    opts.mv_bdries = 0;
end

param_data = chunkerarcparam_init(chnkr);

if ~opts.mv_bdries
    % only shuffle points in each panel
    chnkr = arc_pts(chnkr,param_data);
    eps =  param_data.eps;
else
    % globally reparameterize
    fcurve0 = @(s) chunkerarcparam(s,param_data);
    
    cparams = [];
    cparams.maxchunklen = max(param_data.plen);
    cparams.eps = 10 * param_data.eps;
    
    % get components
    [~,~,info] = sortinfo(chnkr);
    ncomp = info.ncomp;
    nchs = info.nchs;
    
    lens = zeros(1,ncomp);
    for i = 1:ncomp
        lens(i) = sum(chnkr.wts(:,(1:nchs(i)) +sum(nchs(1:i-1))),'all');
    end
    
    comp_len = [0, cumsum(lens)];
    
    % build each component individually, then merge
    chnkrs = [];
    for i = 1:ncomp
        fcurve = @(s) fcurve0(s + comp_len(i));
        cparams.tb = lens(i);
        chnkri = chunkerfunc(fcurve,cparams);
        chnkrs = [chnkrs, chnkri];
    end
    
    chnkr = merge(chnkrs);
    eps = cparams.eps;
end
end


function chnkr = arc_pts(chnkr,param_data)
% parameterize each panel of chunker by arclength

% extract param_data
cr  = param_data.cr;
cd  = param_data.cd;
cd2 = param_data.cd2;
k = param_data.k;
dim = param_data.dim;
nch = param_data.nch;

plen = param_data.plen;

xs = lege.exps(k);
% evaluate legendre series
legs = lege.pols(xs,k-1).';
hs = reshape(plen/2,1,1,[]);
r  = pagetranspose(pagemtimes(legs,cr));
d  = pagetranspose(pagemtimes(legs,cd)).*hs;
d2 = pagetranspose(pagemtimes(legs,cd2)).*hs.^2;

% update chunker
chnkr.rstor(:,:,1:nch) = r;
chnkr.dstor(:,:,1:nch) = d;
chnkr.d2stor(:,:,1:nch) = d2;
chnkr.nstor(:,:,1:nch) = normals(chnkr);
chnkr.wtsstor(:,1:nch) = weights(chnkr);

end
