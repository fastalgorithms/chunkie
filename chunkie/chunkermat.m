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
%           opts.rcip = boolean (true), flag for whether to include rcip
%                      corrections for near corners if input chnkobj is
%                      of type chunkergraph
%           opts.nsub_or_tol = (40) specify the level of refinements in rcip
%                    or a tolerance where the number of levels is given by
%                    ceiling(log_{2}(1/tol^2);
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
isrcip = true;
nsub = 40;

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

if(isfield(opts,'rcip'))
    isrcip = opts.rcip;
end

if(isfield(opts,'nsub_or_tol'))
    if(opts.nsub_or_tol <1)
        tol = opts.nsub_or_tol;
        nsub = max(ceil(log2(1/tol^2)),200);
    else
        nsub = ceil(opts.nsub_or_tol);
    end
    
end

% Flag for determining whether input object is a chunkergraph
icgrph = 0;

if (class(chnkobj) == "chunker")
    chnkrs = chnkobj;
elseif(class(chnkobj) == "chunkgraph")
    icgrph = 1;
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
            ktmp = kern{i,j};
            ftemp = ktmp(srcinfo,targinfo);
        end   
        opdims = size(ftemp);
        opdims_mat(:,i,j) = opdims;
    end
end    

irowlocs = zeros(nchunkers+1,1);
icollocs = zeros(nchunkers+1,1);

irowlocs(1) = 1;
icollocs(1) = 1;
for i=1:nchunkers
   icollocs(i+1) = icollocs(i) + lchunks(i)*opdims_mat(2,1,i);
   irowlocs(i+1) = irowlocs(i) + lchunks(i)*opdims_mat(1,i,1);
end    

nrows = irowlocs(end)-1;
ncols = icollocs(end)-1;

%% Off diagonal interactions 

if (~nonsmoothonly)
    
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

                if (l2scale)
                    wts = sqrt(wts);
                end

                if (size(kern) == 1)
                    ftmp = kern;
                else
                    ftmp = kern{i,j};
                end 

                sysmat_tmp = ftmp(chnkrj,chnkri).*wts;
                
                if (l2scale)
                    wts = weights(chnkri); wts = sqrt(wts(:))';
                    wtsrow = repmat(wts,opdims(1),1); wtsrow = wtsrow(:);
                    sysmat_tmp = bsxfun(@times,wtsrow,sysmat_tmp);
                end
                irowinds = irowlocs(i):(irowlocs(i+1)-1);
                icolinds = icollocs(j):(icollocs(j+1)-1);
                sysmat(irowinds,icolinds) = sysmat_tmp;

           end
        end
    end    
else
    sysmat = sparse(nrows,ncols);
    isysmat = [];
    jsysmat = [];
    vsysmat = [];
end    

%% Diagonal Interaction. 

for i=1:nchunkers
    
    opdims = reshape(opdims_mat(:,i,i),[2,1]);
    jlist = [];
    if ~isempty(ilist)
        jlist = ilist(:,i);
    end

    chnkr = chnkrs(i);
    if (size(kern) == 1)
        ftmp = kern;
    else
        
        ftmp = kern{i,i};
    end 
    
    
    % call requested routine

    if strcmpi(quad,'ggq')
        if (isfield(opts,'auxquads') &&isfield(opts.auxquads,'ggqlog'))
            auxquads = opts.auxquads.ggqlog;
        else
            
        end   
        k = chnkr.k;
        auxquads = chnk.quadggq.setuplogquad(k,opdims);
        opts.auxquads.ggqlog = auxquads;
        type = 'log';
        if nonsmoothonly
            sysmat_tmp = chnk.quadggq.buildmattd(chnkr,ftmp,opdims,type,auxquads,jlist);
        else
            sysmat_tmp = chnk.quadggq.buildmat(chnkr,ftmp,opdims,type,auxquads,jlist);
        end

    elseif strcmpi(quad,'native')

        if nonsmoothonly
            sysmat_tmp = sparse(chnkr.npt,chnkr.npt);
        else
            if (quadorder ~= chnkr.k)
                warning(['native rule: quadorder', ... 
                    ' must equal chunker order (%d)'],chnkr.k)
            end
            sysmat_tmp = chnk.quadnative.buildmat(chnkr,ftmp,opdims);
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
 
    if (~nonsmoothonly)
        irowinds = irowlocs(i):(irowlocs(i+1)-1);
        icolinds = icollocs(i):(icollocs(i+1)-1);
        sysmat(irowinds,icolinds) = sysmat_tmp;
    else
        [isys,jsys,vsys] = find(sysmat_tmp);
        isysmat = [isysmat;isys+irowlocs(i)-1];
        jsysmat = [jsysmat;jsys+icollocs(i)-1];
        vsysmat = [vsysmat;vsys];
    end    
end



if(icgrph && isrcip)
    nch_all = horzcat(chnkobj.echnks.nch);
    [~,nv] = size(chnkobj.verts);
    ngl = chnkrs(1).k;
    glxs = lege.exps(k);
    
    
    for ivert=1:nv
        fprintf('ivert=%d\n',ivert);
        clist = chnkobj.vstruc{ivert}{1};
        isstart = chnkobj.vstruc{ivert}{2};
        isstart(isstart==1) = 0;
        isstart(isstart==-1) = 1;
        nedge = length(isstart);
        iedgechunks = zeros(2,nedge);
        iedgechunks(1,:) = clist;
        iedgechunks(2,:) = 1;
        nch_use = nch_all(clist);
        iedgechunks(2,isstart==0) = nch_use(isstart==0);
        
        
        % since opdims mat for all chunkers meeting at the same 
        % vertex should be the same 
        
        % Todo: have a check that fails with error if this is not the case
        ndim = opdims_mat(1,clist(1),clist(1));
        
        opdims_test = opdims_mat(1:2,clist,clist);
        
        if(norm(opdims_test(:)-ndim)>=0.5) 
            fprintf('in chunkermat: rcip: opdims did not match up for for vertex =%d\n',ivert)
            fprtinf('returning without doing any rcip correction\m');
            break
        end
        starind = zeros(1,2*ngl*ndim*nedge);
        for i=1:nedge
            i1 = (i-1)*2*ngl*ndim+1;
            i2 = i*2*ngl*ndim;
            if(isstart(i))
                starind(i1:i2) = irowlocs(clist(i))+(1:2*ngl*ndim)-1;
            else
                starind(i1:i2) = irowlocs(clist(i)+1)-fliplr(0:2*ngl*ndim-1)-1;
            end
        end
        

        [Pbc,PWbc,starL,circL,starS,circS,ilist] = chnk.rcip.setup(ngl,ndim, ...
          nedge,isstart);
        
        % this might need to be fixed in triple junction case
        R = chnk.rcip.Rcompchunk(chnkrs,iedgechunks,kern,ndim, ...
            Pbc,PWbc,nsub,starL,circL,starS,circS,ilist,... 
            glxs);
       
        sysmat_tmp = inv(R) - eye(2*ngl*nedge*ndim);
        if (~nonsmoothonly)
            
            sysmat(starind,starind) = sysmat_tmp;
        else
            
            isysmat = [isysmat;starind];
            jsysmat = [jsysmat;starind];
            vsysmat = [vsysmat;sysmat_tmp(:)];
        end    
    end
    
end



if (nonsmoothonly)
    sysmat = sparse(isysmat,jsysmat,vsysmat,nrows,ncols);
end


if (nargout >1) 
	varargout{1} = opts;
end  

end
