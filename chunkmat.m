function [sysmat] = chunkmat(chnkr,kern,opts)
%CHUNKMAT build matrix for given kernel and chunker description of 
% boundary. This is a wrapper for various quadrature routines
%
% Input:
%   chnkr - chunker object describing boundary
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(s,t,taus,taut), where s and t
%           are the source and target locations and taus and taut are the 
%           unit tangent at the source and target locations.
%   opts  - options structure. available options (default settings)
%           opts.quad = string ('ggqlog'), specify quadrature routine to 
%                       use. 
%
%                       - 'ggqlog' uses a generalized Gaussian quadrature 
%                       designed for logarithmically singular kernels and 
%                       smooth kernels with removable singularities
%                       - 'native' selects standard scaled Gauss-Legendre 
%                       quadrature for native functions
%
%           opts.quadorder = integer (chnkr.k), desired quadrature order.
%           opts.nonsmoothonly = boolean (false), if true, only compute the
%                         entries for which a special quadrature is used
%                         (e.g. self and neighbor interactoins) and return
%                         in a sparse array.
%           opts.l2scale = boolean (false), if true scale rows by 
%                           sqrt(whts) and columns by 1/sqrt(whts)
%
% Output:
%   sysmat - the system matrix for convolution of the kernel defined by
%            kern with a density on the domain defined by chnkr

if length(chnkr) > 1
    chnkr = chunkermerge(chnkr);
end

if nargin < 3
    opts = [];
end

quadorder = chnkr.k;
quad = 'ggqlog';
nonsmoothonly = false;
l2scale = false;

if or(chnkr.nch < 1,chnkr.k < 1)
    sysmat = [];
    return
end

% determine operator dimensions using first two points

rs = chnkr.r(:,1:2);
ds = chnkr.d(:,1:2); dsn = sqrt(sum(ds.^2,1)); 
ds = bsxfun(@rdivide,ds,dsn);

ftemp = kern(rs(:,1),rs(:,2),ds(:,1),ds(:,2));
opdims = size(ftemp);

% get opts from struct if available

if isfield(opts,'quadorder')
    quadorder = opts.quadorder;
end
if isfield(opts,'quad')
    quad = opts.quad;
end
if isfield(opts,'l2scale')
    l2scale = opts.l2scale;
end
if isfield(opts,'nonsmoothonly')
    nonsmoothonly = opts.nonsmoothonly;
end

% call requested routine

if strcmpi(quad,'ggqlog')
    
    type = 'log';
    if nonsmoothonly
        sysmat = chnk.quadggq.buildmattd(chnkr,kern,opdims,type);
    else
        sysmat = chnk.quadggq.buildmat(chnkr,kern,opdims,type);
    end
    
elseif strcmpi(quad,'native')
        
    if nonsmoothonly
        sysmat = sparse(chnkr.npt,chnkr.npt);
    else
        if (quadorder ~= chnkr.k)
            warning(['native rule: quadorder', ... 
                ' must equal chunker order (%d)'],chnkr.k)
        end
        sysmat = chnk.quadnative.buildmat(chnkr,kern,opdims);
    end
else
    warning('specified quadrature method not available');
    sysmat = [];
    return;
end
	 
if l2scale
    wts = whts(chnkr); wts = sqrt(wts(:)); wts = wts.';
    wtscol = repmat(wts,opdims(2),1); wtscol = wtscol(:); 
    wtscol = wtscol.';
    wtsrow = repmat(wts,opdims(1),1); wtsrow = wtsrow(:);
    sysmat = bsxfun(@times,wtsrow,sysmat);
    sysmat = bsxfun(@rdivide,sysmat,wtscol);
end

end
