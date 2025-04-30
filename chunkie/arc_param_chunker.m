function [chnkr,eps] = arc_param_chunker(chnkr)
% ARC_PARAM_CHUNKER reparameterize chnkr by arclength
%
% NOTE: requires an accurate chunker discretization
% Routine garauntees maxchunklen the does not increase
%
% Output: reparameterized chnkr and tolerance (eps) used to construct it
%
% see also chunkerarcparam

% author: Tristan Goodwill


param_data = chunkerarcparam_init(chnkr);

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