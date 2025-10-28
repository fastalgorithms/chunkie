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
%           opts.eps = (1e-14) tolerance for adaptive quadrature
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
%             for subsequent postprocessing solution at targets closed 
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
            sysmat_tmp = sparse(chnkr.npt,chnkr.npt);
        else
            sysmat_tmp = chnk.quadnative.buildmat(chnkr,ftmp,opdims);
        end
    else
        warning('specified quadrature method not available');
        return;
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

        [R,rcipsav{ivert}] = chnk.rcip.Rcompchunk(chnkrs,iedgechunks,kern,ndim, ...
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
