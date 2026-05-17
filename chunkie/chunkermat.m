function [sysmat,varargout] = chunkermat(chnkobj,kern,opts,ilist)
%CHUNKERMAT build matrix for given kernel and chunker description of 
% boundary. This is a wrapper for various quadrature routines. Optionally,
% return only those interactions which do not use the smooth integration
% rule in the sparse matrix format.
%
% Syntax: sysmat = chunkermat(chnkr,kern,opts)
%
% Input:
%   chnkobj - chunker object describing boundary
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
%                       use for neighbor and self interactoins. Other 
%                       available options include:
%                       - 'native' selects standard scaled Gauss-Legendre 
%                       quadrature for smooth integral kernels
%           opts.sing = string ('log') 
%                       default type of singularity of kernel in case it 
%                       is not defined by the kernel object. Supported 
%                       types are:
%                         smooth => smooth kernels
%                         removable => piecewise smooth kernels
%                         log => logarithmically singular kernels or 
%                                smooth times log + smooth
%                         pv => principal value singular kernels + log
%                         hs => hypersingular kernels + pv
%
%           opts.nonsmoothonly = boolean (false), if true, only compute the
%                         entries for which a special quadrature is used
%                         (e.g. self and neighbor interactions) and return
%                         in a sparse array.
%           opts.corrections = boolean (false), if true, only compute the
%                         corrections to the smooth quadrature rule and 
%                         return in a sparse array, see opts.nonsmoothonly
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
%           opts.rcip_ignore = [], list of vertices to ignore in rcip
%           opts.rcip_savedepth = (10), depth to save rcip info
%           opts.nsub_or_tol = (40) specify the level of refinements in rcip
%                    or a tolerance where the number of levels is given by
%                    ceiling(log_{2}(1/tol^2));
%           opts.adaptive_correction = (false) flag for whether to use
%                    adaptive quadrature for near touching panels on
%                    different chunkers
%           opts.rcip_adaptive_correction = (false) flag for whether to use
%                    adaptive quadrature for near touching panels on
%                    different chunkers within rcip
%           opts.eps = (1e-14) tolerance for adaptive quadrature
%           opts.forcewlchs = kernel-split self/adjacent-panel correction
%                    for Helmholtz / Laplace S, D, S', D' kernels (and
%                    Helmholtz combined 'c' = c1*D + c2*S, 'sc' = c1*S'
%                    + c2*S kernels). After GGQ buildmat fills self+adj
%                    with generic singular quadrature, overwrite
%                    supported entries with wLCHS values from
%                    chnk.kernsplit. Two forms supported:
%                       opts.forcewlchs = true
%                            (kern must be a single helm2d or lap2d
%                            scalar kernel; scalar multiplier auto-
%                            extracted by probing).
%                       opts.forcewlchs = struct('zk', zk, 'layout', L)
%                            where L is an opdims(1)-by-opdims(2) cell
%                            array; L{r,c} = [] for a zero block, or
%                            {type, coef} / struct('type',t,'coef',a)
%                            for a scalar kernel. zk required for
%                            Helmholtz family entries.
%  ilist - cell array of integer arrays ([]), list of panel interactions that
%          should be ignored when constructing matrix entries or quadrature
%          corrections. 
%
% Output:
%   sysmat - the system matrix for discretizing integral operator whose kernel 
%            is defined by kern with a density on the domain defined by chnkr
%
% Optional output
%   opts - with the updated opts structure which stores the relevant
%          quantities in opts.auxquads.<opts.quad><opts.type>
%   rcipsav - precomputed structure of rcip data at corners
%             for subsequent postprocessing of the solution at targets close 
%             to the corner
%
% Examples:
%   sysmat = chunkermat(chnkr,kern); % standard options
%   sysmat = chunkermat(chnkr,kern,opts);
%   sysmat = chunkermat(chnkr,kern,opts,ilist);
%   [sysmat,opts] = chunkermat(chnkr,kern,opts);
%   [sysmat,opts] = chunkermat(chnkr,kern,opts,ilist);
%

% convert kernel to kernel object, put in singularity info 
% opts.sing provides a default value for singularities if not 
% defined for kernels


% Flag for determining whether input object is a chunkergraph
icgrph = 0;

if (class(chnkobj) == "chunker")
    chnkrs = chnkobj;
    npttot = chnkobj.npt;
elseif(class(chnkobj) == "chunkgraph")
    icgrph = 1;
    chnkrs = chnkobj.echnks;
    npttot = chnkobj.npt;
else
    msg = "CHUNKERMAT: first input is not a chunker or chunkgraph object";
    error(msg)
end

if ~isa(kern,'kernel')
    try 
        kern = kernel(kern);
    catch
        error('CHUNKERMAT: second input kern not of supported type');
    end
end

if nargin < 3
    opts = [];
end

if nargin <4
    ilist = [];
end

quad = 'ggq';
nonsmoothonly = false;
corrections = false;
l2scale = false;
isrcip = true;
rcip_ignore = [];
nsub = 40;
rcip_savedepth = 10;
adaptive_correction = false;
rcip_adaptive_correction = false;
sing = 'log';

% get opts from struct if available

if isfield(opts,'quad')
    quad = opts.quad;
end
if isfield(opts,'sing')
    sing = opts.sing;
end
if isfield(opts,'l2scale')
    l2scale = opts.l2scale;
end
if isfield(opts,'nonsmoothonly')
    nonsmoothonly = opts.nonsmoothonly;
end
if isfield(opts,'corrections')
    corrections = opts.corrections;
end
if corrections
    nonsmoothonly = true;
end

if(isfield(opts,'rcip'))
    isrcip = opts.rcip;
end
if(isfield(opts,'rcip_savedepth'))
    rcip_savedepth = opts.rcip_savedepth;
end
if(isfield(opts,'rcip_ignore'))
    rcip_ignore = opts.rcip_ignore;
    if ~isrcip && isempty(rcip_ignore)
        fprintf('in chunkermat: provided list of vertices to ignore in RCIP\n')
        fprintf('while RCIP is not enabled.\n')
    end
end

if(isfield(opts,'nsub_or_tol'))
    if(opts.nsub_or_tol <1)
        tol = opts.nsub_or_tol;
        nsub = max(ceil(log2(1/tol^2)),200);
    else
        nsub = ceil(opts.nsub_or_tol);
    end
end

if (isfield(opts,'adaptive_correction'))
    adaptive_correction = opts.adaptive_correction;
end
if (isfield(opts,'rcip_adaptive_correction'))
    rcip_adaptive_correction = opts.rcip_adaptive_correction;
end

% forcewlchs: kernel-split self/adjacent-panel correction for Laplace
% S/D/Sp/Dp kernels. See doc above for accepted forms.
forcewlchs_opt = [];
if isfield(opts,'forcewlchs') && ~isempty(opts.forcewlchs)
    forcewlchs_opt = opts.forcewlchs;
end

nchunkers = length(chnkrs);

opdims_mat = zeros(2,nchunkers,nchunkers);
lchunks    = zeros(nchunkers,1);

%TODO: figure out a way to avoid this nchunkers^2 loop

for i=1:nchunkers
    
    targinfo = [];
   	targinfo.r = chnkrs(i).r(:,2); targinfo.d = chnkrs(i).d(:,2); 
   	targinfo.d2 = chnkrs(i).d2(:,2); targinfo.n = chnkrs(i).n(:,2);
    if ~isempty(chnkrs(i).data)
	    targinfo.data = chnkrs(i).data(:, 2);
    end
    lchunks(i) = chnkrs(i).npt;
    
    for j=1:nchunkers
        
        % determine operator dimensions using first two points

        srcinfo = []; 
        srcinfo.r = chnkrs(j).r(:,1); srcinfo.d = chnkrs(j).d(:,1); 
        srcinfo.d2 = chnkrs(j).d2(:,1); srcinfo.n = chnkrs(j).n(:,1);
        if ~isempty(chnkrs(j).data)
	        srcinfo.data = chnkrs(j).data(:, 1);
        end

        if (size(kern) == 1)
            ftemp = kern.eval(srcinfo,targinfo);
        else
            ktmp = kern(i,j).eval;
            ftemp = ktmp(srcinfo,targinfo);
        end   
        opdims = size(ftemp);
        opdims_mat(:,i,j) = opdims;
    end
end    

irowlocs = zeros(nchunkers+1,1);
icollocs = zeros(nchunkers+1,1);

idrowchnk = zeros(2,npttot);
idcolchnk = zeros(2,npttot);


irowlocs(1) = 1;
icollocs(1) = 1;
for i=1:nchunkers
   icollocs(i+1) = icollocs(i) + lchunks(i)*opdims_mat(2,1,i);
   irowlocs(i+1) = irowlocs(i) + lchunks(i)*opdims_mat(1,i,1);

    % which chunker
    idrowchnk(1,irowlocs(i):(irowlocs(i+1)-1)) = i;
    % which chunk in the chunker
    idrowchnk(2,irowlocs(i):(irowlocs(i+1)-1)) = ceil((1:lchunks(i)*opdims_mat(1,i,1))/ ...
        (opdims_mat(1,i,1)*chnkrs(i).k));

    % which chunker
    idcolchnk(1,icollocs(i):(icollocs(i+1)-1)) = i;
    % which chunk in the chunker
    idcolchnk(2,icollocs(i):(icollocs(i+1)-1)) = ceil((1:lchunks(i)*opdims_mat(2,1,i))/ ...
        (opdims_mat(2,1,i)*chnkrs(i).k));
end    

nrows = irowlocs(end)-1;
ncols = icollocs(end)-1;

if (~nonsmoothonly)
    sysmat = zeros(nrows,ncols);
else 
    sysmat = sparse(nrows,ncols);
    isysmat = [];
    jsysmat = [];
    vsysmat = [];
end


%% Off diagonal interactions 

%TODO: switch to flagging for all chunks at once, then loop over
% chunks and do corrections. need to create array from chunk/point index
% to edge chunk number

for i = 1:nchunkers
    chnkri = chnkrs(i);
    for j = 1:nchunkers
        chnkrj = chnkrs(j);
        if (chnkri.nch < 1 || chnkri.k < 1 || chnkrj.nch<1 || chnkrj.k<1)
            sysmat_tmp = [];
            break
        end

        if (i~=j)

            opdims = reshape(opdims_mat(:,i,j),[2,1]);
            wts = chnkrj.wts;
            wts2 = repmat( (wts(:)).', opdims(2), 1);
            wts2 = ( wts2(:) ).';
            wts = wts2;

            if (size(kern) == 1)
                ftmp = kern.eval;
            else
                ftmp = kern(i,j).eval;
            end 

            if (~nonsmoothonly)
                sysmat_tmp = ftmp(chnkrj,chnkri).*wts;
            end

            % Cross-chunker (i ~= j) near-corner correction via wLCHS.
            % When forcewlchs is enabled and the (i,j) block kernel is a
            % recognized scalar/combined Helmholtz type, add the kernel-
            % split close-evaluation delta from helm2d_near_correction.
            % The proximity check inside that function filters out non-
            % near (i,j) pairs, so chunker pairs that don't meet at a
            % corner contribute zero.
            if ~isempty(forcewlchs_opt) && ~nonsmoothonly && all(opdims == 1)
                if size(kern,1) == 1 && size(kern,2) == 1
                    bk = kern;
                else
                    bk = kern(i,j);
                end
                [near_type, near_coef] = chunkermat_kernel_for_near(bk);
                if ~isempty(near_type)
                    if any(strcmp(near_type, {'c','combined','sc','spcombined'}))
                        delta_near = chnk.kernsplit.helm2d_near_correction( ...
                            chnkrj, chnkri, near_type, bk.params.zk, [], near_coef);
                    else
                        delta_near = chnk.kernsplit.helm2d_near_correction( ...
                            chnkrj, chnkri, near_type, bk.params.zk);
                        delta_near = near_coef * delta_near;
                    end
                    if nnz(delta_near) > 0
                        sysmat_tmp = sysmat_tmp + full(delta_near);
                    end
                end
            end

            if adaptive_correction
                flag = flagnear(chnkrj,chnkri.r(:,:));
                sysmat_tmp_adap = chunkermat_adap(chnkrj,ftmp,opdims, ...
                    chnkri,flag,opts,corrections);
                if (~nonsmoothonly)
                    [isys,jsys,vsys] = find(sysmat_tmp_adap);
                    ijsys = isys + (jsys-1)*size(sysmat_tmp_adap,1);
                    sysmat_tmp(ijsys) = vsys(:);
                else
                    sysmat_tmp = sysmat_tmp_adap;

                end
            end

            if l2scale
                wts = sqrt(wts); 
                wtsrow = chnkri.wts; wtsrow = sqrt(wtsrow(:))';
                wtsrow = repmat(wtsrow,opdims(1),1); wtsrow = wtsrow(:);
                sysmat_tmp = wtsrow.*sysmat_tmp./wts;
            end
            
            if (~nonsmoothonly)
                irowinds = irowlocs(i):(irowlocs(i+1)-1);
                icolinds = icollocs(j):(icollocs(j+1)-1);
                sysmat(irowinds,icolinds) = sysmat_tmp;
            else
                if adaptive_correction
                    [isys,jsys,vsys] = find(sysmat_tmp);
                    isysmat = [isysmat;isys+irowlocs(i)-1];
                    jsysmat = [jsysmat;jsys+icollocs(j)-1];
                    vsysmat = [vsysmat;vsys];
                end
            end
        end
    end
end

%% Diagonal Interaction. 

for i=1:nchunkers
    
    opdims = reshape(opdims_mat(:,i,i),[2,1]);
    jlist = [];
    if ~isempty(ilist)
        jlist = ilist(:,i);
    end

    singi = sing;
    chnkr = chnkrs(i);
    if (size(kern) == 1)
        ftmp = kern.eval;
        if ~isempty(kern.sing)
            singi = kern.sing;
        end
    else
        
        ftmp = kern(i,i).eval;
        if ~isempty(kern(i,i).sing)
            singi = kern(i,i).sing;
        end
    end 

    
    
    % call requested routine

    if strcmpi(quad,'ggq')
        if strcmpi(singi,'smooth')
            %TODO: make a reasonable method for smooth with removable
            type = 'log';
            if (isfield(opts,'auxquads') &&isfield(opts.auxquads,'ggqlog'))
                auxquads = opts.auxquads.ggqlog;
            else
                k = chnkr.k;
                auxquads = chnk.quadggq.setup(k,type);
                opts.auxquads.ggqlog = auxquads;
            end
        elseif strcmpi(sing, 'removable')
            type = 'removable';
            if (isfield(opts,'auxquads') && isfield(opts.auxquads,'ggqremovable'))
                auxquads = opts.auxquads.ggqremovable;
            else
                k = chnkr.k;
                auxquads = chnk.quadggq.setup(k, type);
                opts.auxquads.ggqremovable = auxquads;
            end
        elseif strcmpi(singi,'log')
            type = 'log';
            if (isfield(opts,'auxquads') &&isfield(opts.auxquads,'ggqlog'))
                auxquads = opts.auxquads.ggqlog;
            else
                k = chnkr.k;
                auxquads = chnk.quadggq.setup(k,type);
                opts.auxquads.ggqlog = auxquads;
            end
        elseif strcmpi(singi,'pv')
            type = 'pv';
            if (isfield(opts,'auxquads') &&isfield(opts.auxquads,'ggqpv'))
                auxquads = opts.auxquads.ggqpv;
            else
                k = chnkr.k;
                auxquads = chnk.quadggq.setup(k,type);
                opts.auxquads.ggqpv = auxquads;
            end
        elseif strcmpi(singi,'hs')
            type = 'hs';
            if (isfield(opts,'auxquads') &&isfield(opts.auxquads,'ggqhs'))
                auxquads = opts.auxquads.ggqhs;
            else
                k = chnkr.k;
                auxquads = chnk.quadggq.setup(k,type);
                opts.auxquads.ggqhs = auxquads;
            end
        end
        if nonsmoothonly
            sysmat_tmp = chnk.quadggq.buildmattd(chnkr,ftmp,opdims,type,auxquads,jlist,corrections);
        else
            sysmat_tmp = chnk.quadggq.buildmat(chnkr,ftmp,opdims,type,auxquads,jlist);
        end

    elseif strcmpi(quad,'native')

        if nonsmoothonly
            sysmat_tmp = sparse(opdims(1)*chnkr.npt, opdims(2)*chnkr.npt);
        else
            sysmat_tmp = chnk.quadnative.buildmat(chnkr,ftmp,opdims);
        end
    else
        warning('specified quadrature method not available');
        return;
    end

    if ~isempty(forcewlchs_opt)
        % In ~nonsmoothonly mode, sysmat_tmp has smooth GL values for the
        % whole block (native quad) or smooth + GGQ self/adj corrections
        % (ggq quad). chunkermat_apply_wlchs_adj ADDs the kernel-split
        % delta over smooth GL on self+adj entries. When ggq quad is in
        % use, the GGQ self/adj corrections it wrote get overwritten by
        % the wLCHS values -- so prefer quad='native' for cleanliness.
        % In nonsmoothonly mode, sysmat_tmp is sparse (zero from native,
        % or GGQ corrections from ggq); add wLCHS deltas as sparse too.
        if size(kern,1) > 1 || size(kern,2) > 1
            block_kern = kern(i,i);
        else
            block_kern = kern;
        end
        try
            sysmat_tmp = chunkermat_apply_wlchs_adj( ...
                sysmat_tmp, chnkr, opdims, block_kern, forcewlchs_opt, nonsmoothonly);
        catch ME
            if contains(ME.message, 'forcewlchs') ...
                    || contains(ME.message, 'wlchs') ...
                    || contains(ME.message, 'unsupported')
                warning('chunkermat:wlchsSkipped', ...
                    'forcewlchs skipped on chunker %d: %s', i, ME.message);
            else
                rethrow(ME);
            end
        end
    end

    if adaptive_correction
        flag = flagnear(chnkr,chnkr.r(:,:));

        % mark off the near and self interactions
        for ich = 1:chnkr.nch
	    jlist = [ich,chnkr.adj(1,ich),chnkr.adj(2,ich)];
	    jlist = jlist(jlist > 0);
            for jch = jlist
                flag((jch - 1)*chnkr.k+(1:chnkr.k), ich) = 0;
            end
        end
        
        sysmat_tmp_adap = chunkermat_adap(chnkr, ftmp, opdims, chnkr, ...
           flag,opts,corrections);

        [isys,jsys,vsys] = find(sysmat_tmp_adap);

        ijsys = isys + (jsys-1)*size(sysmat_tmp_adap,1);
        sysmat_tmp(ijsys) = vsys(:);
    end

    if l2scale
        wts = chnkr.wts; wts = sqrt(wts(:)); wts = wts.';
        wtscol = repmat(wts,opdims(2),1); wtscol = wtscol(:); 
        wtscol = wtscol.';
        wtsrow = repmat(wts,opdims(1),1); wtsrow = wtsrow(:);
        sysmat_tmp = wtsrow.*sysmat_tmp./wtscol;
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
    [sbclmat,sbcrmat,lvmat,rvmat,u] = chnk.rcip.shiftedlegbasismats(k);
    nch_all = horzcat(chnkobj.echnks.nch);
    npt_all = horzcat(chnkobj.echnks.npt);
    [~,nv] = size(chnkobj.verts);
    ngl = chnkrs(1).k;

    rcipsav = cell(nv,1);
    
    for ivert=setdiff(1:nv,rcip_ignore)
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
            fprintf('returning without doing any rcip correction\n');
            break
        end
        starind = zeros(1,2*ngl*ndim*nedge);
        corinds = cell(nedge,1);
        for i=1:nedge
            i1 = (i-1)*2*ngl*ndim+1;
            i2 = i*2*ngl*ndim;
            if(isstart(i))
                starind(i1:i2) = irowlocs(clist(i))+(1:2*ngl*ndim)-1;
                corinds{i} = 1:2*ngl;
            else
                starind(i1:i2) = irowlocs(clist(i)+1)-fliplr(0:2*ngl*ndim-1)-1;
                corinds{i} = (npt_all(clist(i))-2*ngl + 1):npt_all(clist(i));
            end
        end
        
        [Pbc,PWbc,starL,circL,starS,circS,ilist,starL1,circL1] = ...
            chnk.rcip.setup(ngl,ndim,nedge,isstart);
        optsrcip = opts;
        optsrcip.nonsmoothonly = false;
        optsrcip.corrections = false;
        optsrcip.rcip_savedepth = rcip_savedepth;
        optsrcip.adaptive_correction = rcip_adaptive_correction;


        [R,rcipsav{ivert}] = chnk.rcip.Rcompchunk(chnkrs,iedgechunks,kern,ndim,chnkobj.verts(:,ivert), ...
            Pbc,PWbc,nsub,starL,circL,starS,circS,ilist,starL1,circL1,... 
            sbclmat,sbcrmat,lvmat,rvmat,u,optsrcip);

        rcipsav{ivert}.starind = starind;

        sysmat_tmp = inv(R) - eye(2*ngl*nedge*ndim);
        if (~nonsmoothonly)
            
            sysmat(starind,starind) = sysmat_tmp;
        else

            if corrections
                cormat = zeros(size(sysmat_tmp));
                jstart = 1;
                for jj = 1:nedge
                    jch = clist(jj);
                    srcinfo = [];
                    jinds = corinds{jj};
                    srcinfo.r = chnkrs(jch).r(:,jinds);
                    srcinfo.d = chnkrs(jch).d(:,jinds);
                    srcinfo.d2 = chnkrs(jch).d2(:,jinds);
                    srcinfo.n = chnkrs(jch).n(:,jinds);
                    wtsj = chnkrs(jch).wts(jinds);
                    op2 = opdims_mat(2,1,jch);
                    wtsj = repmat(wtsj(:).',op2,1);
                    wtsj = wtsj(:);
                    jflat = jstart:(jstart-1+length(jinds)*op2);
                    jstart = jstart+length(jinds)*op2;
                    istart = 1;
                    for ii = 1:nedge
                        ich = clist(ii);
                        iinds = corinds{ii};
                        targinfo = [];
                        targinfo.r = chnkrs(ich).r(:,iinds);
                        targinfo.d = chnkrs(ich).d(:,iinds);
                        targinfo.d2 = chnkrs(ich).d2(:,iinds);
                        targinfo.n = chnkrs(ich).n(:,iinds);
                        op1 = opdims_mat(1,ich,1);
                        iflat = istart:(istart-1+length(iinds)*op1);
                        istart = istart + length(iinds)*op1;
                        if size(kern) == 1
                            ktmp = kern.eval;
                        else
                            ktmp = kern(ich,jch).eval;
                        end
                        submat = ktmp(srcinfo,targinfo).*(wtsj(:).');
                        if (ii == jj)
                            submat(kron(eye(length(iinds)),ones(op1,op2)) > 0) = 0;
                        end
                        cormat(iflat,jflat) = submat;
                    end
                end
                sysmat_tmp = sysmat_tmp - cormat;
            end
                       
            [jind,iind] = meshgrid(starind);
            
            isysmat = [isysmat;iind(:)];
            jsysmat = [jsysmat;jind(:)];
            vsysmat = [vsysmat;sysmat_tmp(:)];
        end    
    end

    if nargout > 2
        varargout{2} = rcipsav;
    end

    
end




if (nonsmoothonly)
    % Fix sparse entry format to use rcip matrix entries for repeats
    % instead of using the precomputed self correction
    
    ijind = [isysmat jsysmat];
    [~,idx] = unique(ijind, 'rows','last');
    
    sysmat = sparse(isysmat(idx),jsysmat(idx),vsysmat(idx),nrows,ncols);
end


if (nargout >1) 
	varargout{1} = opts;
end  

end

function mat = chunkermat_adap(chnkr,kern,opdims, ...
			       chnkrt,flag,opts,corrections)

  if isa(kern,'kernel')
    kernev = kernel.eval;
  else
    kernev = kern;
  end

% copied and modified from the function chunkerkernevalmat in chunkerkernevalmat.m
% chnkrt is the target chunker. 

k = chnkr.k;
nch = chnkr.nch;

if nargin < 5
    flag = [];
end
if nargin < 6
    opts = [];
end
if nargin < 7
    corrections = false;
end

targs = chnkrt.r(:,:); targn = chnkrt.n(:,:); 
targd = chnkrt.d(:,:); targd2 = chnkrt.d2(:,:);
targdata = chnkrt.data(:,:);

[~,nt] = size(targs);

if (~any(flag(:)))
    mat = sparse(opdims(1)*nt,opdims(2)*chnkr.npt);
    return
end

% using adaptive quadrature

is = zeros(nnz(flag)*opdims(1)*opdims(2)*k,1);
js = is;
vs = is;
istart = 1;

[t,w] = lege.exps(2*k+1);
ct = chnkr.tstor;
bw = lege.barywts(k,ct);
r = chnkr.r;
d = chnkr.d;
n = chnkr.n;
d2 = chnkr.d2;
data = chnkr.data;

wtss = chnkr.wts;

for i = 1:nch
    jmat = 1 + (i-1)*k*opdims(2);
    jmatend = i*k*opdims(2);
                    
    [ji] = find(flag(:,i));
    if ~isempty(targdata)
        mat1 =  chnk.adapgausswts(r,d,n,d2,data,ct,bw,i,targs(:,ji), ...
                targd(:,ji),targn(:,ji),targd2(:,ji),targdata(:,ji),...
                kernev,opdims,t,w,opts);
    else
        mat1 =  chnk.adapgausswts(r,d,n,d2,data,ct,bw,i,targs(:,ji), ...
                targd(:,ji),targn(:,ji),targd2(:,ji),[],...
                kernev,opdims,t,w,opts);
    end            
    if corrections
        targinfo = []; targinfo.r = targs(:,ji); targinfo.d = targd(:,ji);
        targinfo.n = targn(:,ji); targinfo.d2 = targd2(:,ji);
        if ~isempty(targdata)
            targinfo.data = targdata(:,ji);
        end
        srcinfo = []; srcinfo.r = r(:,:,i); srcinfo.d = d(:,:,i); 
        srcinfo.n = n(:,:,i); srcinfo.d2 = d2(:,:,i);
        if ~isempty(data)
            srcinfo.data = data(:,:,i);
        end
        wtsi = wtss(:,i); wtsi = repmat(wtsi(:).',opdims(2),1);
        mat1 = mat1 - kernev(srcinfo,targinfo).*(wtsi(:).');
    end
    js1 = jmat:jmatend;
    js1 = repmat( (js1(:)).',opdims(1)*numel(ji),1);
            
    indji = (ji-1)*opdims(1);
    indji = repmat( (indji(:)).', opdims(1),1) + ( (1:opdims(1)).');
    indji = indji(:);
    
    indji = repmat(indji,1,opdims(2)*k);
    
    iend = istart+numel(mat1)-1;
    is(istart:iend) = indji(:);
    js(istart:iend) = js1(:);
    vs(istart:iend) = mat1(:);
    istart = iend+1;
end

mat = sparse(is,js,vs,opdims(1)*nt,opdims(2)*chnkr.npt);

end


function sysmat = chunkermat_apply_wlchs_adj(sysmat, chnkr, opdims, kern, fw, nonsmoothonly)
% Apply wLCHS self/adjacent-panel corrections to sysmat for Helmholtz /
% Laplace S/D/Sp/Dp entries. See chunkermat opts.forcewlchs comment.
%
% Behaviour depends on nonsmoothonly:
%   false (default) - sysmat is dense; ENTRIES on self+adj blocks are
%                     OVERWRITTEN with smooth-GL + wLCHS-delta. Used by
%                     the dense matrix-build path.
%   true            - sysmat is sparse; wLCHS deltas (smooth-subtracted)
%                     are ADDED on top of whatever GGQ corrections were
%                     already in sysmat. For matrix-free apply, the caller
%                     should also use quad='native' so sysmat starts at
%                     zero and the result is just the wLCHS delta.
if nargin < 6 || isempty(nonsmoothonly), nonsmoothonly = false; end

% Decode forcewlchs spec: build a 2D cell `layout` (opdims(1) x opdims(2))
% where layout{r,c} is either [] (no wLCHS, leave entry alone), or
% struct('type', t, 'coef', a, 'zk', zk, 'family', 'helmholtz'|'laplace').
% Helmholtz types: s, d, sp, dp, t, c, sc. Laplace types: s, d, sp, dp
% (no combined kernels). zk is unused/ignored for Laplace.
m = opdims(1); n = opdims(2);
layout = cell(m, n);
helm_types = {'s','single','d','double','sp','sprime','dp','dprime', ...
              'c','combined','sc','spcombined'};
lap_types  = {'s','single','d','double','sp','sprime','dp','dprime'};
if islogical(fw) || (isnumeric(fw) && isscalar(fw))
    % Boolean form: kern must be a single scalar Helmholtz or Laplace
    % kernel, or a zero kernel (no-op).
    if ~fw, return; end
    if isa(kern, 'kernel') && (strcmpi(kern.name, 'zeros') || ...
            (isprop(kern, 'iszero') && kern.iszero))
        return;
    end
    fam = wlchs_kernel_family(kern);
    if isempty(fam)
        error("chunkermat: opts.forcewlchs=true requires a single helm2d or lap2d kernel");
    end
    t = lower(kern.type);
    if m ~= 1 || n ~= 1
        error("chunkermat: opts.forcewlchs=true requires a single scalar (opdims=[1 1]) helm or lap kernel");
    end
    if strcmp(fam, 'helmholtz')
        if ~ismember(t, helm_types)
            error("chunkermat: opts.forcewlchs=true: helm kernel type %s not supported", kern.type);
        end
        if ~isfield(kern.params, 'zk') || isempty(kern.params.zk)
            error("chunkermat: opts.forcewlchs requires kern.params.zk");
        end
        zk = kern.params.zk;
    else  % laplace
        if ~ismember(t, lap_types)
            error("chunkermat: opts.forcewlchs=true: lap kernel type %s not supported", kern.type);
        end
        zk = [];
    end
    if strcmpi(t, 'c') || strcmpi(t, 'combined')
        if ~isfield(kern.params, 'coefs') || isempty(kern.params.coefs)
            error("chunkermat: forcewlchs c kernel requires kern.params.coefs");
        end
        layout{1,1} = struct('type', 'c', 'coefs', kern.params.coefs(:).', 'zk', zk, 'family', fam);
    elseif strcmpi(t, 'sc') || strcmpi(t, 'spcombined')
        if ~isfield(kern.params, 'coefs') || isempty(kern.params.coefs)
            error("chunkermat: forcewlchs sc kernel requires kern.params.coefs");
        end
        layout{1,1} = struct('type', 'sc', 'coefs', kern.params.coefs(:).', 'zk', zk, 'family', fam);
    else
        coef = wlchs_extract_scalar(kern, t, zk, fam);
        layout{1,1} = struct('type', t, 'coef', coef, 'zk', zk, 'family', fam);
    end
elseif isstruct(fw)
    if ~isfield(fw, 'layout')
        error("chunkermat: opts.forcewlchs struct needs field 'layout'");
    end
    % Top-level defaults; per-entry struct can override family / zk.
    fam_default = 'helmholtz';
    if isfield(fw, 'family') && ~isempty(fw.family)
        fam_default = lower(fw.family);
    end
    if ~ismember(fam_default, {'helmholtz','laplace'})
        error("chunkermat: opts.forcewlchs.family must be 'helmholtz' or 'laplace'");
    end
    zk_default = [];
    if isfield(fw, 'zk') && ~isempty(fw.zk), zk_default = fw.zk; end

    L = fw.layout;
    if ~iscell(L) || size(L,1) ~= m || size(L,2) ~= n
        error("chunkermat: opts.forcewlchs.layout must be %dx%d cell", m, n);
    end
    for r = 1:m
        for c = 1:n
            entry = L{r,c};
            if isempty(entry); continue; end
            ent_fam = fam_default;
            ent_zk  = zk_default;
            if iscell(entry)
                t = entry{1}; a = entry{2};
            elseif isstruct(entry)
                t = entry.type;
                if isfield(entry, 'coefs') && ~isempty(entry.coefs)
                    a = entry.coefs;
                else
                    a = entry.coef;
                end
                if isfield(entry, 'family') && ~isempty(entry.family)
                    ent_fam = lower(entry.family);
                end
                if isfield(entry, 'zk') && ~isempty(entry.zk), ent_zk = entry.zk; end
            else
                error("chunkermat: opts.forcewlchs.layout{%d,%d} bad form", r, c);
            end
            if strcmpi(t, 'zero'); continue; end

            switch ent_fam
                case 'helmholtz'
                    allowed_e = {'s','single','d','double','sp','sprime','dp','dprime', ...
                                 't','c','combined','sc','spcombined'};
                    if isempty(ent_zk)
                        error("chunkermat: forcewlchs helm layout entry (%d,%d) needs zk", r, c);
                    end
                case 'laplace'
                    allowed_e = lap_types;
                otherwise
                    error("chunkermat: unknown family %s", ent_fam);
            end
            if ~ismember(lower(t), allowed_e)
                error("chunkermat: forcewlchs (%s): layout entry type %s not supported", ent_fam, t);
            end

            if strcmpi(t, 'c') || strcmpi(t, 'combined')
                if ~isvector(a) || numel(a) ~= 2
                    error("chunkermat: forcewlchs c layout entry needs coefs [c1, c2]");
                end
                layout{r,c} = struct('type', 'c', 'coefs', a(:).', ...
                                     'zk', ent_zk, 'family', ent_fam);
            elseif strcmpi(t, 'sc') || strcmpi(t, 'spcombined')
                if ~isvector(a) || numel(a) ~= 2
                    error("chunkermat: forcewlchs sc layout entry needs coefs [c1, c2]");
                end
                layout{r,c} = struct('type', 'sc', 'coefs', a(:).', ...
                                     'zk', ent_zk, 'family', ent_fam);
            else
                layout{r,c} = struct('type', lower(t), 'coef', a, ...
                                     'zk', ent_zk, 'family', ent_fam);
            end
        end
    end
else
    error("chunkermat: opts.forcewlchs has unsupported type");
end

% For each non-empty (r,c) in layout, OVERWRITE self and adjacent block
% entries with smooth-GL + wLCHS-delta (replacing the inadequate generic
% GGQ values that the buildmat already wrote).
ngl = chnkr.k;
% U_src (per-panel inverse Vandermonde for wlchs) is kernel-INDEPENDENT,
% so build it once across the layout iteration. Cells fill lazily as
% lap2d_*_correction touch them.
U_src_cell = cell(1, chnkr.nch);
for r = 1:m
    for c = 1:n
        spec = layout{r,c};
        if isempty(spec); continue; end

        is_c  = strcmpi(spec.type, 'c')  || strcmpi(spec.type, 'combined');
        is_sc = strcmpi(spec.type, 'sc') || strcmpi(spec.type, 'spcombined');
        is_split = is_c || is_sc;

        % --- Self-panel correction ---
        % Plain scalar kernels: spec.coef * delta_self(spec.type)
        % Combined 'c' = c1*D + c2*S:   c2 * delta_self('s') for helm
        %                               (Helm D has no self correction).
        % Combined 'sc' = c1*S' + c2*S: c1*delta_self('sp') + c2*delta_self('s')
        % (Combined kernels are Helmholtz-only.)
        if wlchs_has_self_correction(spec.type, spec.family)
            if is_c
                c1 = spec.coefs(1); c2 = spec.coefs(2);
                delta_self_a = []; delta_self_b = [];
                if c2 ~= 0
                    [delta_self_a, U_src_cell] = wlchs_self_correction_dispatch( ...
                        spec, chnkr, 's', U_src_cell);
                end
                self_coef_a = c2; self_coef_b = 0;
            elseif is_sc
                c1 = spec.coefs(1); c2 = spec.coefs(2);
                delta_self_a = []; delta_self_b = [];
                if c1 ~= 0
                    [delta_self_a, U_src_cell] = wlchs_self_correction_dispatch( ...
                        spec, chnkr, 'sp', U_src_cell);
                end
                if c2 ~= 0
                    [delta_self_b, U_src_cell] = wlchs_self_correction_dispatch( ...
                        spec, chnkr, 's', U_src_cell);
                end
                self_coef_a = c1; self_coef_b = c2;
            else
                [delta_self_a, U_src_cell] = wlchs_self_correction_dispatch( ...
                    spec, chnkr, spec.type, U_src_cell);
                delta_self_b = [];
                self_coef_a = spec.coef; self_coef_b = 0;
            end
            for ich = 1:chnkr.nch
                ji = (1:ngl).' + (ich-1)*ngl;
                rows_blk = (ji - 1) * m + r;
                cols_blk = (ji.' - 1) * n + c;
                d_blk = zeros(ngl, ngl);
                if ~isempty(delta_self_a)
                    d_blk = d_blk + self_coef_a * full(delta_self_a(ji, ji));
                end
                if ~isempty(delta_self_b)
                    d_blk = d_blk + self_coef_b * full(delta_self_b(ji, ji));
                end
                if nonsmoothonly
                    sysmat(rows_blk, cols_blk) = d_blk;
                else
                    if is_split
                        Kb = wlchs_kernel_eval(spec, spec.type, ...
                            chnkr.r(:,:,ich), chnkr.n(:,:,ich), ...
                            chnkr.r(:,:,ich), chnkr.n(:,:,ich), spec.coefs);
                        awzp = chnkr.wts(:, ich).';
                        M_smooth = (Kb .* awzp);
                    else
                        Kb = wlchs_kernel_eval(spec, spec.type, ...
                            chnkr.r(:,:,ich), chnkr.n(:,:,ich), ...
                            chnkr.r(:,:,ich), chnkr.n(:,:,ich));
                        awzp = chnkr.wts(:, ich).';
                        M_smooth = spec.coef * (Kb .* awzp);
                    end
                    M_smooth(1:ngl+1:end) = 0;
                    sysmat(rows_blk, cols_blk) = M_smooth + d_blk;
                end
            end
        end

        % --- Adjacent-panel correction ---
        delta_a = []; delta_b = [];
        adj_coef_a = 0; adj_coef_b = 0;
        if is_c
            c1 = spec.coefs(1); c2 = spec.coefs(2);
            if c1 ~= 0
                [delta_a, U_src_cell] = wlchs_adj_correction_dispatch( ...
                    spec, chnkr, 'd', U_src_cell);
            end
            if c2 ~= 0
                [delta_b, U_src_cell] = wlchs_adj_correction_dispatch( ...
                    spec, chnkr, 's', U_src_cell);
            end
            adj_coef_a = c1; adj_coef_b = c2;
        elseif is_sc
            c1 = spec.coefs(1); c2 = spec.coefs(2);
            if c1 ~= 0
                [delta_a, U_src_cell] = wlchs_adj_correction_dispatch( ...
                    spec, chnkr, 'sp', U_src_cell);
            end
            if c2 ~= 0
                [delta_b, U_src_cell] = wlchs_adj_correction_dispatch( ...
                    spec, chnkr, 's', U_src_cell);
            end
            adj_coef_a = c1; adj_coef_b = c2;
        else
            [delta_a, U_src_cell] = wlchs_adj_correction_dispatch( ...
                spec, chnkr, spec.type, U_src_cell);
            adj_coef_a = spec.coef;
        end

        for ich = 1:chnkr.nch
            for io = 1:2
                jch = chnkr.adj(io, ich);
                if jch <= 0; continue; end
                ji    = (1:ngl).' + (ich-1)*ngl;
                jcols = (jch-1)*ngl + (1:ngl);
                rows_blk = (ji - 1) * m + r;
                cols_blk = (jcols - 1) * n + c;
                d_blk = zeros(ngl, ngl);
                if ~isempty(delta_a)
                    d_blk = d_blk + adj_coef_a * full(delta_a(ji, jcols));
                end
                if ~isempty(delta_b)
                    d_blk = d_blk + adj_coef_b * full(delta_b(ji, jcols));
                end
                if nonsmoothonly
                    sysmat(rows_blk, cols_blk) = d_blk;
                else
                    if is_split
                        Kb = wlchs_kernel_eval(spec, spec.type, ...
                            chnkr.r(:,:,jch), chnkr.n(:,:,jch), ...
                            chnkr.r(:,:,ich), chnkr.n(:,:,ich), spec.coefs);
                        awzp = chnkr.wts(:, jch).';
                        M_smooth = (Kb .* awzp);
                    else
                        Kb = wlchs_kernel_eval(spec, spec.type, ...
                            chnkr.r(:,:,jch), chnkr.n(:,:,jch), ...
                            chnkr.r(:,:,ich), chnkr.n(:,:,ich));
                        awzp = chnkr.wts(:, jch).';
                        M_smooth = spec.coef * (Kb .* awzp);
                    end
                    sysmat(rows_blk, cols_blk) = M_smooth + d_blk;
                end
            end
        end
    end
end
end

function fam = wlchs_kernel_family(kern)
% Returns 'helmholtz', 'laplace', or '' if the kernel is unsupported by
% the wLCHS path.
fam = '';
if ~isa(kern, 'kernel'), return; end
if strcmpi(kern.name, 'zeros'), return; end
if strcmpi(kern.name, 'helmholtz'), fam = 'helmholtz'; return; end
if strcmpi(kern.name, 'laplace'),   fam = 'laplace';   return; end
end

function tf = wlchs_has_self_correction(type, family)
% Helmholtz: S, S', D', T, combined kernels have self corrections; D does
% not (the J_1 piece of D vanishes on the diagonal at smooth boundary
% points). For combined 'c' = c1*D + c2*S we return true and let the
% application step multiply c2 (zero if S is absent); same idea for 'sc'.
% Laplace: S, D, S', D' all have nontrivial self corrections.
if nargin < 2 || isempty(family), family = 'helmholtz'; end
switch lower(family)
    case 'helmholtz'
        tf = ismember(lower(type), {'s','single','sp','sprime','dp','dprime','t', ...
                                    'c','combined','sc','spcombined'});
    case 'laplace'
        tf = ismember(lower(type), {'s','single','d','double','sp','sprime','dp','dprime'});
    otherwise
        tf = false;
end
end

function [delta, U_src_cell] = wlchs_self_correction_dispatch(spec, chnkr, type, U_src_cell)
% Dispatch self-panel correction to the right kernel family. The `type`
% argument is passed explicitly because combined kernels (Helm 'c'/'sc')
% issue per-component calls with types differing from spec.type.
switch lower(spec.family)
    case 'helmholtz'
        [delta, U_src_cell] = chnk.kernsplit.helm2d_self_correction(chnkr, type, spec.zk, U_src_cell);
    case 'laplace'
        [delta, U_src_cell] = chnk.kernsplit.lap2d_self_correction(chnkr, type, U_src_cell);
    otherwise
        error('wlchs_self_correction_dispatch: unsupported family %s', spec.family);
end
end

function [delta, U_src_cell] = wlchs_adj_correction_dispatch(spec, chnkr, type, U_src_cell)
switch lower(spec.family)
    case 'helmholtz'
        [delta, U_src_cell] = chnk.kernsplit.helm2d_adj_correction(chnkr, type, spec.zk, U_src_cell);
    case 'laplace'
        [delta, U_src_cell] = chnk.kernsplit.lap2d_adj_correction(chnkr, type, U_src_cell);
    otherwise
        error('wlchs_adj_correction_dispatch: unsupported family %s', spec.family);
end
end

function K = wlchs_kernel_eval(spec, type, src_r, src_n, tgt_r, tgt_n, varargin)
% Evaluate the bare kernel of given family/type at (src, tgt) GL nodes.
% For combined types, varargin{1} = coefs = [c1, c2] (Helmholtz only).
si.r = src_r; si.n = src_n;
ti.r = tgt_r; ti.n = tgt_n;
switch lower(spec.family)
    case 'helmholtz'
        zk = spec.zk;
        switch lower(type)
            case {'s','single'},      K = chnk.helm2d.kern(zk, si, ti, 's');
            case {'d','double'},      K = chnk.helm2d.kern(zk, si, ti, 'd');
            case {'sp','sprime'},     K = chnk.helm2d.kern(zk, si, ti, 'sprime');
            case {'dp','dprime','t'}, K = chnk.helm2d.kern(zk, si, ti, 'dprime');
            case {'c','combined'}
                if isempty(varargin)
                    error('wlchs_kernel_eval: c kernel requires coefs');
                end
                K = chnk.helm2d.kern(zk, si, ti, 'c', varargin{1});
            case {'sc','spcombined'}
                if isempty(varargin)
                    error('wlchs_kernel_eval: sc kernel requires coefs');
                end
                K = chnk.helm2d.kern(zk, si, ti, 'sc', varargin{1});
            otherwise
                error('wlchs_kernel_eval: helm type %s not supported', type);
        end
    case 'laplace'
        switch lower(type)
            case {'s','single'},  K = chnk.lap2d.kern(si, ti, 's');
            case {'d','double'},  K = chnk.lap2d.kern(si, ti, 'd');
            case {'sp','sprime'}, K = chnk.lap2d.kern(si, ti, 'sprime');
            case {'dp','dprime'}, K = chnk.lap2d.kern(si, ti, 'dprime');
            otherwise
                error('wlchs_kernel_eval: lap type %s not supported', type);
        end
    otherwise
        error('wlchs_kernel_eval: unsupported family %s', spec.family);
end
end

function coef = wlchs_extract_scalar(kern, t, param, family)
% Probe the kern eval at a known src/tgt pair, divide by analytic value.
% Result is the scalar prefactor on top of the bare kernel (=1 for plain
% kernel('lap',type) or kernel('helm',type,zk); !=1 for scalar-multiplied
% kernels). For Laplace S the bare value at r=1 is 0 (log 1 = 0), so we
% probe at r=2 instead.
%
% `param` carries family-specific data (zk for helmholtz, [] for laplace).
if nargin < 4 || isempty(family), family = 'helmholtz'; end
if strcmp(family, 'laplace')
    src_probe = struct('r',[0;0],'n',[1;0],'d',[1;0],'d2',[0;0]);
    tgt_probe = struct('r',[2;0],'n',[1;0],'d',[1;0],'d2',[0;0]);
else
    src_probe = struct('r',[0;0],'n',[1;0],'d',[1;0],'d2',[0;0]);
    tgt_probe = struct('r',[1;0],'n',[1;0],'d',[1;0],'d2',[0;0]);
end
val = kern.eval(src_probe, tgt_probe);
switch lower(family)
    case 'helmholtz'
        zk = param;
        switch lower(t)
            case {'s','single'}
                base = 0.25i * besselh(0, zk);
            case {'d','double'}
                base = 0.25i * zk * besselh(1, zk);
            case {'sp','sprime'}
                base = kernel('helm','sp', zk).eval(src_probe, tgt_probe);
            case {'dp','dprime'}
                base = kernel('helm','dp', zk).eval(src_probe, tgt_probe);
            otherwise
                error('wlchs_extract_scalar: helm type %s not supported', t);
        end
        coef = val / base;
    case 'laplace'
        switch lower(t)
            case {'s','single'},  base = chnk.lap2d.kern(src_probe, tgt_probe, 's');
            case {'d','double'},  base = chnk.lap2d.kern(src_probe, tgt_probe, 'd');
            case {'sp','sprime'}, base = chnk.lap2d.kern(src_probe, tgt_probe, 'sprime');
            case {'dp','dprime'}, base = chnk.lap2d.kern(src_probe, tgt_probe, 'dprime');
            otherwise
                error('wlchs_extract_scalar: lap type %s not supported', t);
        end
        coef = val / base;
    otherwise
        error('wlchs_extract_scalar: unsupported family %s', family);
end
end


function [near_type, near_coef] = chunkermat_kernel_for_near(kern)
% Inspect a kernel object and return (type, coef-or-coefs) suitable for
% calling chnk.kernsplit.helm2d_near_correction. Returns empty type for
% unsupported kernels (zero kernel, custom types, etc.).
near_type = ''; near_coef = [];
if ~isa(kern, 'kernel'); return; end
if strcmpi(kern.name, 'zeros') || (isprop(kern,'iszero') && kern.iszero)
    return;
end
if ~strcmpi(kern.name, 'helmholtz'); return; end
if ~isfield(kern.params, 'zk') || isempty(kern.params.zk); return; end
zk_l = kern.params.zk;
t = lower(kern.type);
switch t
    case {'s','single','d','double','sp','sprime','dp','dprime','t'}
        near_type = t;
        near_coef = wlchs_extract_scalar(kern, t, zk_l, 'helmholtz');
    case {'c','combined'}
        if isfield(kern.params, 'coefs') && numel(kern.params.coefs) == 2
            near_type = 'c';
            near_coef = kern.params.coefs(:).';
        end
    otherwise
        % unsupported (custom sum kernels with no recognized type)
end
end
