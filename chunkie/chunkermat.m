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

if nargin < 3
    opts = [];
end

if nargin <4
    ilist = [];
end

quad = 'ggq';
nonsmoothonly = false;
l2scale = false;

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


if (class(chnkobj) == "chunker")
    chnkrs = chnkobj;
elseif(class(chnkobj) == "chunkgraph")
    chnkrs = chnkobj.echnks;
else
    msg = "Unsupported object in chunkermat";
    error(msg)
end

nchunkers = length(chnkrs);

opdims_mat = zeros(2,nchunkers,nchunkers);
lchunks    = zeros(nchunkers,1);

for i=1:nchunkers
    
    targinfo = [];
   	targinfo.r = chnkrs(i).r(:,2); targinfo.d = chnkrs(i).d(:,2); 
   	targinfo.d2 = chnkrs(i).d2(:,2); targinfo.n = chnkrs(i).n(:,2);
    lchunks(i) = size(chnkrs(i).r(:,:),2);
    
    for j=1:nchunkers
        
        % determine operator dimensions using first two points

        srcinfo = []; 
        srcinfo.r = chnkrs(j).r(:,1); srcinfo.d = chnkrs(j).d(:,1); 
        srcinfo.d2 = chnkrs(j).d2(:,1); srcinfo.n = chnkrs(j).n(:,1);

        if (size(kern) == 1)
            ftemp = kern(srcinfo,targinfo);
        else
            ktmp = kern(i,j);
            ftemp = ktmp(srcinfo,targinfo);
        end   
        opdims = size(ftemp);
        opdims_mat(:,i,j) = opdims;
    end
end    

irowlocs = [1];
icollocs = [1];

for i=1:nchunkers
   icollocs(i+1) = icollocs(i) + lchunks(i)*opdims_mat(2,1,i);
   irowlocs(i+1) = irowlocs(i) + lchunks(i)*opdims_mat(1,i,1);
end    

nrows = irowlocs(end)-1;
ncols = icollocs(end)-1;

sysmat = zeros(nrows,ncols);

for i = 1:nchunkers
    chnkri = chnkrs(i);
    for j = 1:nchunkers
        chnkrj = chnkrs(j);
        if (chnkri.nch < 1 || chnkri.k < 1 || chnkrj.nch<1 || chnkri.k<1)
            sysmat_tmp = [];
            break
        end
        
       if (i~=j)

            opdims = reshape(opdims_mat(:,i,j),[2,1]);
            wts = weights(chnkrj);
            wts2 = repmat( (wts(:)).', opdims(2), 1);
            wts2 = ( wts2(:) ).';
            wts = wts2;

            if (size(kern) == 1)
                ftmp = kern;
            else
                ftmp = kern(i,j);
            end 
            sysmat_tmp = ftmp(chnkrj,chnkri).*wts;
            irowinds = irowlocs(i):(irowlocs(i+1)-1);
            icolinds = icollocs(j):(icollocs(j+1)-1);
            sysmat(irowinds,icolinds) = sysmat_tmp;

       end
    end
end    

for i=1:nchunkers

    opdims = reshape(opdims_mat(:,i,i),[2,1]);

    chnkr = chnkrs(i);
    
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
            sysmat_tmp = chnk.quadggq.buildmattd(chnkr,kern,opdims,type,auxquads,ilist);
        else
            sysmat_tmp = chnk.quadggq.buildmat(chnkr,kern,opdims,type,auxquads,ilist);
        end

    elseif strcmpi(quad,'native')

        if nonsmoothonly
            sysmat_tmp = sparse(chnkr.npt,chnkr.npt);
        else
            if (quadorder ~= chnkr.k)
                warning(['native rule: quadorder', ... 
                    ' must equal chunker order (%d)'],chnkr.k)
            end
            sysmat_tmp = chnk.quadnative.buildmat(chnkr,kern,opdims);
        end
    else
        warning('specified quadrature method not available');
        sysmat_tmp = [];
        return;
    end

    if l2scale
        wts = weights(chnkr); wts = sqrt(wts(:)); wts = wts.';
        wtscol = repmat(wts,opdims(2),1); wtscol = wtscol(:); 
        wtscol = wtscol.';
        wtsrow = repmat(wts,opdims(1),1); wtsrow = wtsrow(:);
        sysmat_tmp = bsxfun(@times,wtsrow,sysmat_tmp);
        sysmat_tmp = bsxfun(@rdivide,sysmat_tmp,wtscol);
    end
 
    irowinds = irowlocs(i):(irowlocs(i+1)-1);
  	icolinds = icollocs(i):(icollocs(i+1)-1);
	sysmat(irowinds,icolinds) = sysmat_tmp;

end



if (nargout >1) 
	varargout{1} = opts;
end  

end
