function [sysmat,varargout] = chunkermat(chnkobj,kern,opts,ilist)
%CHUNKERMAT build matrix for given kernel and chunker description of 
% boundary. This is a wrapper for various quadrature routines. Optionally,
% return only those interactions which do not use the smooth integration
% rule in the sparse matrix format.
%
% Syntax: sysmat = chunkermat(chnkr,kern,opts)
%
% Input:
%   chnkr - chunker object describing boundary
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(srcinfo,targinfo), where srcinfo
%           and targinfo are in the ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.n - unit normals (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%
% Optional input:
%   opts  - options structure. available options (default settings)
%           opts.quad = string ('ggq'), specify quadrature routine to 
%                       use. Other available options include
%                       - 'native' selects standard scaled Gauss-Legendre 
%                       quadrature for native functions
%                       smooth kernels with removable singularities
%           opts.type = string ('log'), type of singularity of kernel. Type
%                       can take on the following arguments:
%                         log => logarithmically singular kernels
%                         pv => principal value singular kernels
%                         hs => hypersingular kernels
%
%           opts.nonsmoothonly = boolean (false), if true, only compute the
%                         entries for which a special quadrature is used
%                         (e.g. self and neighbor interactoins) and return
%                         in a sparse array.
%           opts.l2scale = boolean (false), if true scale rows by 
%                           sqrt(whts) and columns by 1/sqrt(whts)
%           opts.auxquads = struct, struct storing auxilliary nodes 
%                     and weights which might be required for some of
%                     the quadrature methods like ggq for example.
%                     There is a different sub structure for each
%                     quadrature and singularity type which should be named
%                     as
%
%                     opts.auxquads.<opts.quad><opts.type> 
%                     
%                     For example, the structure for logarithmically
%                     singular kernels integrated using ggq
%                     quadrature, the relevant struct is
%                     
%                     opts.auxquads.ggqlog
%
%                     The specific precomputed variables and their values
%                     will depend on the quadrature method used.
%  ilist - cell array of integer arrays ([]), list of panel interactions that 
%          should be ignored when constructing matrix entries or quadrature
%          corrections. 
%
%
% Output:
%   sysmat - the system matrix for discretizing integral operator whose kernel 
%            is defined by kern with a density on the domain defined by chnkr
%
% Optional output
%   opts - with the updated opts structure which stores the relevant
%          quantities in opts.auxquads.<opts.quad><opts.type>
%
% Examples:
%   sysmat = chunkermat(chnkr,kern); % standard options
%   sysmat = chunkermat(chnkr,kern,opts);
%   sysmat = chunkermat(chnkr,kern,opts,ilist);
%   [sysmat,opts] = chunkermat(chnkr,kern,opts);
%   [sysmat,opts] = chunkermat(chnkr,kern,opts,ilist);
%

if (class(chnkobj) == "chunker")
    chnkr = chnkobj;
elseif(class(chnkobj) == "chunkgraph")
    chnkr = merge(chnkobj.echnks);
else
    msg = "Unsupported object in chunkermat";
    error(msg)
end

if length(chnkr) > 1
    chnkr = merge(chnkr);
end

if nargin < 3
    opts = [];
end

if nargin <4
    ilist = [];
end

quad = 'ggq';
nonsmoothonly = false;
l2scale = false;

if or(chnkr.nch < 1,chnkr.k < 1)
    sysmat = [];
    return
end

% determine operator dimensions using first two points

srcinfo = []; targinfo = [];
srcinfo.r = chnkr.r(:,1); srcinfo.d = chnkr.d(:,1); 
srcinfo.d2 = chnkr.d2(:,1); srcinfo.n = chnkr.n(:,1);
targinfo.r = chnkr.r(:,2); targinfo.d = chnkr.d(:,2); 
targinfo.d2 = chnkr.d2(:,2); targinfo.n = chnkr.n(:,2);

ftemp = kern(srcinfo,targinfo);
opdims = size(ftemp);

% get opts from struct if available

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

if strcmpi(quad,'ggq')
    if (isfield(opts,'auxquads') &&isfield(opts.auxquads,'ggqlog'))
        auxquads = opts.auxquads.ggqlog;
    else
        k = chnkr.k;
        auxquads = chnk.quadggq.setuplogquad(k,opdims);
        opts.auxquads.ggqlog = auxquads;
    end    
    type = 'log';
    if nonsmoothonly
        sysmat = chnk.quadggq.buildmattd(chnkr,kern,opdims,type,auxquads,ilist);
    else
        sysmat = chnk.quadggq.buildmat(chnkr,kern,opdims,type,auxquads,ilist);
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
    wts = weights(chnkr); wts = sqrt(wts(:)); wts = wts.';
    wtscol = repmat(wts,opdims(2),1); wtscol = wtscol(:); 
    wtscol = wtscol.';
    wtsrow = repmat(wts,opdims(1),1); wtsrow = wtsrow(:);
    sysmat = bsxfun(@times,wtsrow,sysmat);
    sysmat = bsxfun(@rdivide,sysmat,wtscol);
end

if (nargout >1) 
   varargout{1} = opts;
end   

end
