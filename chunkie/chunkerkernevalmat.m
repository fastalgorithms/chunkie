function [mat,varargout] = chunkerkernevalmat(chnkobj,kern,targobj,opts)
%CHUNKERKERNEVALMAT compute the matrix which maps density values on 
% the chunk geometry to the value of the convolution of the given
% integral kernel with the density at the specified target points
%
% Syntax: mat = chunkerkernevalmat(chnkobj,kern,targs,opts)
%
% Input:
%   chnkobj - chunker object or chunkgraph object description of curve
%   kern  - kernel function. By default, this should be a function handle
%           accepting input of the form kern(srcinfo,targinfo), where srcinfo
%           and targinfo are in the ptinfo struct format, i.e.
%                ptinfo.r - positions (2,:) array
%                ptinfo.d - first derivative in underlying
%                     parameterization (2,:)
%                ptinfo.n - unit normals (2,:)
%                ptinfo.d2 - second derivative in underlying
%                     parameterization (2,:)
%   targobj - object describing the target points, can be specified as
%       * array of points
%       * chunker object
%       * chunkgraph object
%
% Optional input:
%   opts - structure for setting various parameters
%       opts.forcesmooth - if = true, only use the smooth integration rule
%                           (false)
%       opts.forceadap - if = true, only use adaptive quadrature (false)
%    NOTE: only one of forcesmooth or forceadap is allowed. If both 
%           false, a hybrid algorithm is used, where the smooth rule is 
%           applied for targets separated by opts.fac*length of chunk from
%           a given chunk and adaptive integration is used otherwise
%       opts.fac = the factor times the chunk length used to decide 
%               between adaptive/smooth rule
%       opts.eps = tolerance for adaptive integration
%       opts.nonsmoothonly = boolean (false), if true, only compute the
%                         entries for which a special quadrature is used
%                         (e.g. self and neighbor interactions) and return
%                         in a sparse array.
%       opts.corrections = boolean (false), if true, only compute the
%                         corrections to the smooth quadrature rule and 
%                         return in a sparse array, see opts.nonsmoothonly

%
% output:
%   mat - (opdims(1)*nt) x opdims(2)*chnkr.npt matrix mapping values 
%         of a density on the boundary to values of convolution at target 
%         points
%
% see also CHUNKERKERNEVAL


% author: Travis Askham (askhamwhat@gmail.com)

% Assign appropriate object to chnkr
icgrph = false;
nregion = 1;
nedge = 1;
if class(chnkobj) == "chunker"
   chnkr = chnkobj;
   nedge = length(chnkr);
elseif class(chnkobj) == "chunkgraph"
   chnkr = chnkobj.echnks;
   nregion = length(chnkobj.regions);
   nedge = length(chnkr);
   icgrph = true;
else
    msg = "CHUNKERKERNEVAL: first input is an unsupported object";
    error(msg)
end

if ~isa(kern,'kernel')
    try 
        kern = kernel(kern);
    catch
        error('CHUNKERKERNEVALMAT: second input kern not of supported type');
    end
end

[mk,nk] = size(kern);
assert(or(mk == 1,mk == nregion),...
    "CHUNKERKERNEVAL: second input not of appropriate shape " + ...
    "number of rows in kern should be 1 or nregion")
assert(or(nk == 1,nk == nedge),...
    "CHUNKERKERNEVAL: second input not of appropriate shape " + ...
    "number of cols in kern should be 1 or nedge")

if nk == 1 && length(chnkr) > 1
    chnkr = merge(chnkr);
end
    
if nargin < 4
    opts = [];
end

forcesmooth = false;
forceadap = false;
forcepquad = false;
nonsmoothonly = false;
corrections = false;
fac = 1.0;
eps = 1e-12;
if isfield(opts,'forcesmooth'); forcesmooth = opts.forcesmooth; end
if isfield(opts,'forceadap'); forceadap = opts.forceadap; end
if isfield(opts,'forcepquad'); forcepquad = opts.forcepquad; end
if isfield(opts,'nonsmoothonly'); nonsmoothonly = opts.nonsmoothonly; end
if isfield(opts,'corrections'); corrections = opts.corrections; end
if corrections; nonsmoothonly = true; forcesmooth=false; end

if isfield(opts,'fac'); fac = opts.fac; end
if isfield(opts,'eps'); eps = opts.eps; end

% Assign appropriate object to targinfo
targinfo = [];
if isa(targobj, "chunker")
    targinfo.r = targobj.r(:,:);
    targinfo.d = targobj.d(:,:);
    targinfo.d2 = targobj.d2(:,:);
    targinfo.n = targobj.n(:,:);
    targinfo.data = targobj.data(:,:);
elseif isa(targobj, "chunkgraph")
    targinfo.r = targobj.r(:,:);
    targinfo.d = targobj.d(:,:);
    targinfo.d2 = targobj.d2(:,:);
    targinfo.n = targobj.n(:,:);
    targinfo.data = targinfo.data(:,:);
elseif isstruct(targobj)
    if isfield(targobj,"r")
        targinfo.r = targobj.r(:,:);
    else
        error("CHUNKERKERNEVAL: input 4 must at least have positions " + ...
            "defined");
    end
    if isfield(targobj,"d"); targinfo.d = targobj.d(:,:); end
    if isfield(targobj,"d2"); targinfo.d2 = targobj.d2(:,:); end
    if isfield(targobj,"n"); targinfo.n = targobj.n(:,:); end
    if isfield(targobj,"data"); targinfo.data = targobj.data(:,:); end
elseif isnumeric(targobj)
    targinfo.r = targobj;
else
    error("CHUNKERKERNEVAL: input 4 is not a supported type");
end

if icgrph && mk > 1
    ids = chunkgraphinregion(chnkobj,targinfo.r);
else
    ids = ones(size(targinfo.r,2),1);
end

[dim,~] = size(targinfo.r);

if (dim ~= 2); warning('only dimension two tested'); end

opdims_mat = zeros(2,mk,nk);
ntargs = zeros(mk,1);
npts = zeros(nk,1);

datadim = 0;
if isfield(targinfo,'data') && ~isempty(targinfo.data)
    datadim = size(targinfo.data,1);
end

for iii=1:mk
    itarg = (ids == iii);    
    ntargs(iii) = nnz(itarg);

    targinfotmp = [];
    targinfotmp.r = randn(dim,1); targinfotmp.d = randn(dim,1);
    targinfotmp.d2 = randn(dim,1); targinfotmp.n = randn(dim,1);
    targinfotmp.data = randn(datadim,1);
    
    for jjj=1:nk
        
        % determine operator dimensions using a boundary point and random
        % targ
        
        srcinfo = []; 
        srcinfo.r = chnkr(jjj).r(:,1); srcinfo.d = chnkr(jjj).d(:,1); 
        srcinfo.d2 = chnkr(jjj).d2(:,1); srcinfo.n = chnkr(jjj).n(:,1);
        if ~isempty(chnkr(jjj).data)
            srcinfo.data = chnkr(jjj).data(:,1);
        end
        npts(jjj) = chnkr(jjj).npt; 

        try
            ftemp = kern(iii,jjj).eval(srcinfo,targinfotmp);
        catch
            error("failed to determine size of kernel (%d, %d)",iii,jjj);
        end
        opdims = size(ftemp);
        opdims_mat(:,iii,jjj) = opdims;
    end
end    

% indexing

icollocs = zeros(nk+1,1);
icollocs(1)=1;
for jjj=1:nk
    icollocs(jjj+1) = icollocs(jjj) + npts(jjj)*opdims_mat(2,1,jjj);
end

rowdims = opdims_mat(1,:,1); rowdims = rowdims(:);
nout = sum(ntargs(:).*rowdims(:));

ntarg = size(targinfo.r(:,:),2);
itargstart = zeros(ntarg+1,1);
itargstart(2:end) = rowdims(ids(:));
itargstart = 1+cumsum(itargstart);

optssmooth = []; 
optsadap = []; 
optsadap.eps = eps;

if corrections
    mat = sparse(nout,icollocs(end)-1);
else
    mat = zeros(nout,icollocs(end)-1);
end
for iii = 1:mk
% loop over relevant regions 
itarg = (ids == iii);
if nnz(itarg) == 0
    continue
end

targinfo0 = [];
targinfo0.r = targinfo.r(:,itarg);
if isfield(targinfo,"d"); targinfo0.d = targinfo.d(:,itarg); end
if isfield(targinfo,"d2"); targinfo0.d2 = targinfo.d2(:,itarg); end
if isfield(targinfo,"n"); targinfo0.n = targinfo.n(:,itarg); end
if isfield(targinfo,"data") && ~isempty(targinfo.data); targinfo0.data = targinfo.data(:,itarg); end

irow0 = kron(itargstart(itarg),ones(rowdims(iii),1)) + repmat( (0:(rowdims(iii)-1)).',nnz(itarg),1);

for jjj = 1:nk
% loop over relevant boundary components
icol0 = icollocs(jjj):(icollocs(jjj+1)-1);

kern0 = kern(iii,jjj);
if kern0.isnan
    mat(irow0,icol0) = nan;
    continue
end
if kern0.iszero
    continue
end

chnkr0 = chnkr(jjj);
opdims0 = opdims_mat(:,iii,jjj);

if forcesmooth
    mat(irow0,icol0) = chunkerkernevalmat_smooth(chnkr0,kern0,opdims0,targinfo0, ...
        [],optssmooth);
    continue
end




if forceadap
    mat(irow0,icol0) = chunkerkernevalmat_adap(chnkr0,kern0,opdims0, ...
        targinfo0,[],optsadap);
    continue
end

optsflag = []; optsflag.fac = fac;
flag = flagnear(chnkr0,targinfo0.r,optsflag);

if forcepquad
    spmat = chunkerkernevalmat_pquad(chnkr0,kern0,opdims0, ...
        targinfo0,flag,opts);
    if corrections
        mat(irow0,icol0) = spmat;
        continue
    else
        mat(irow0,icol0) = chunkerkernevalmat_smooth(chnkr0,kern0,opdims0, ...
            targinfo0,flag,opts);
        mat(irow0,icol0) = mat(irow0,icol0) + spmat;
        continue
    end
end

% smooth for sufficiently far, adaptive otherwise

%rho = 1.8;
%optsflag = [];  optsflag.rho = rho;
%flag = flagnear_rectangle(chnkr0,targinfo0.r,optsflag);

%npoly = chnkr.k*2;
%nlegnew = chnk.ellipse_oversample(rho,npoly,eps);
%nlegnew = max(nlegnew,chnkr.k);


spmat = chunkerkernevalmat_adap(chnkr0,kern0,opdims0, ...
        targinfo0,flag,optsadap);


if corrections
    % TODO: find more elegant solution that avoids building a dense flag matrix
    flaginv = ~flag;
    mat0 = spmat - chunkerkernevalmat_smooth2(chnkr0,kern0,opdims0,targinfo0, ...
        flag,opts);
    mat(irow0,icol0) = sparse(mat0);
    continue
end

if nonsmoothonly
    mat(irow0,icol0) = spmat;
    continue;
end

mat(irow0,icol0) = chunkerkernevalmat_smooth(chnkr0,kern0,opdims0,targinfo0, ...
    flag,opts);

mat(irow0,icol0) = mat(irow0,icol0) + spmat;

end
end

end



function mat = chunkerkernevalmat_smooth(chnkr,kern,opdims, ...
    targinfo,flag,opts)

if isa(kern,'kernel')
    kerneval = kern.eval;
else
    kerneval = kern;
end

if nargin < 6
    flag = [];
end
if nargin < 7
    opts = [];
end

k = chnkr.k;
nch = chnkr.nch;

srcinfo = []; srcinfo.r = chnkr.r(:,:); srcinfo.n = chnkr.n(:,:);
srcinfo.d = chnkr.d(:,:); srcinfo.d2 = chnkr.d2(:,:);

mat = kerneval(srcinfo,targinfo);
wts = chnkr.wts;
wts2 = repmat( (wts(:)).', opdims(2), 1);
wts2 = ( wts2(:) ).';
mat = mat.*wts2;

if isempty(flag)
    % nothing to erase
    return
else
    % delete interactions in flag array
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);

        rowkill = find(flag(:,i)); 
        rowkill = (opdims(1)*(rowkill(:)-1)).' + (1:opdims(1)).';
        mat(rowkill,jmat:jmatend) = 0;
    end
end

end

function mat = chunkerkernevalmat_smooth2(chnkr,kern,opdims, ...
    targinfo,flag,opts)

if isa(kern,'kernel')
    kerneval = kern.eval;
else
    kerneval = kern;
end

k = chnkr.k;
nch = chnkr.nch;

if nargin < 5
    flag = [];
end
if nargin < 6
    opts = [];
end

% Extract target info
targs = targinfo.r;
[~,nt] = size(targs);
targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
targn = zeros(chnkr.dim,nt);
if isfield(targinfo, 'd')
    targd = targinfo.d;
end

if isfield(targinfo, 'd2')
    targd2 = targinfo.d2;
end

if isfield(targinfo, 'n')
    targn = targinfo.n;
end
datat = [];
if isfield(targinfo, 'data')
    datat = targinfo.data;
end


% using adaptive quadrature


if isempty(flag)


    mat = sparse(opdims(1)*nt,opdims(2)*chnkr.npt);
    
else
    is = zeros(nnz(flag)*opdims(1)*opdims(2)*k,1);
    js = is;
    vs = is;
    istart = 1;

    r = chnkr.r;
    d = chnkr.d;
    n = chnkr.n;
    d2 = chnkr.d2;
    data = chnkr.data;
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        [ji] = find(flag(:,i));
        datat2 = [];
        if ~isempty(datat)
            datat2 = datat(:,ji);
        end
        srcinfo = [];
        srcinfo.r = r(:,:,i);
        srcinfo.d = d(:,:,i);
        srcinfo.n = n(:,:,i);
        srcinfo.d2 = d2(:,:,i);

        if ~isempty(data); srcinfo.data = data(:,:,i); end

        targinfo = [];
        targinfo.r = targs(:,ji);
        targinfo.d = targd(:,ji);
        targinfo.n = targn(:,ji);
        targinfo.d2 = targd2(:,ji);
        targinfo.data = datat2;

        mat1 = kerneval(srcinfo,targinfo);
        wts1 = chnkr.wts(:,i);
        wts2 = repmat( (wts1(:)).', opdims(2), 1);
        wts2 = ( wts2(:) ).';
        mat1 = mat1.*wts2;

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

end

function mat = chunkerkernevalmat_adap(chnkr,kern,opdims, ...
    targinfo,flag,opts)

if isa(kern,'kernel')
    kerneval = kern.eval;
else
    kerneval = kern;
end

k = chnkr.k;
nch = chnkr.nch;

if nargin < 5
    flag = [];
end
if nargin < 6
    opts = [];
end

% Extract target info
targs = targinfo.r;
[~,nt] = size(targs);
targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
targn = zeros(chnkr.dim,nt);
if isfield(targinfo, 'd')
    targd = targinfo.d;
end

if isfield(targinfo, 'd2')
    targd2 = targinfo.d2;
end

if isfield(targinfo, 'n')
    targn = targinfo.n;
end
datat = [];
if isfield(targinfo, 'data')
    datat = targinfo.data;
end


% using adaptive quadrature


if isempty(flag)


    mat = zeros(opdims(1)*nt,opdims(2)*chnkr.npt);

    [t,w] = lege.exps(2*k+1);
    ct = lege.exps(k);
    bw = lege.barywts(k);
    r = chnkr.r;
    d = chnkr.d;
    n = chnkr.n;
    d2 = chnkr.d2;
    data = chnkr.data;
    
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        mat(:,jmat:jmatend) =  chnk.adapgausswts(r,d,n,d2,data,ct,bw,i,targs, ...
                    targd,targn,targd2,datat,kerneval,opdims,t,w,opts);
    end
    
else
    is = zeros(nnz(flag)*opdims(1)*opdims(2)*k,1);
    js = is;
    vs = is;
    istart = 1;

    [t,w] = lege.exps(2*k+1);
    ct = lege.exps(k);
    bw = lege.barywts(k);
    r = chnkr.r;
    d = chnkr.d;
    n = chnkr.n;
    d2 = chnkr.d2;
    data = chnkr.data;
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        [ji] = find(flag(:,i));
        datat2 = [];
        if ~isempty(datat)
            datat2 = datat(:,ji);
        end
        mat1 =  chnk.adapgausswts(r,d,n,d2,data,ct,bw,i,targs(:,ji), ...
                    targd(:,ji),targn(:,ji),targd2(:,ji),datat2,kerneval,opdims,t,w,opts);
                
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

end

function mat = chunkerkernevalmat_pquad(chnkr,kern,opdims, ...
    targinfo,flag,opts)

if isa(kern,'kernel')
    kerneval = kern.eval;
else
    kerneval = kern;
end

k = chnkr.k;
nch = chnkr.nch;

if nargin < 5
    flag = [];
end
if nargin < 6
    opts = [];
end

% Extract target info
targs = targinfo.r;
[~,nt] = size(targs);
targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
targn = zeros(chnkr.dim,nt);
if isfield(targinfo, 'd')
    targd = targinfo.d;
end

if isfield(targinfo, 'd2')
    targd2 = targinfo.d2;
end

if isfield(targinfo, 'n')
    targn = targinfo.n;
end

dclosest = Inf;

% using Helsing-Ojala quadrature
if isempty(flag) % figure out what is this flag for in adaptive routine
    keyboard
else
    is = zeros(nnz(flag)*opdims(1)*opdims(2)*k,1);
    js = is;
    vs = is;
    istart = 1;

    [t,w] = lege.exps(2*k);
    ct = lege.exps(k);
    bw = lege.barywts(k);
    r = chnkr.r;
    d = chnkr.d;
    n = chnkr.n;
    d2 = chnkr.d2;
    wts = chnkr.wts;

    % interpolation matrix 
    intp = lege.matrin(k,t);          % interpolation from k to 2*k
    intp_ab = lege.matrin(k,[-1;1]);  % interpolation from k to end points
    
    optsuse = opts;
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        [ji] = find(flag(:,i));

        targinfoji = [];
        targinfoji.r = targinfo.r(:,ji);

        srcinfo = [];
        srcinfo.r = r(:,:,i);
        srcinfo.d = d(:,:,i);
        srcinfo.d2 = d2(:,:,i);
        srcinfo.n = n(:,:,i);

        if isfield(opts,'side')
            if strcmp(opts.side,'i')
                iiin = true(1,size(targinfoji.r,2));
                iout = false(1,size(targinfoji.r,2));
            else
                iiin = false(1,size(targinfoji.r,2));
                iout = true(1,size(targinfoji.r,2));
            end
        else
            jjclose = zeros(1,size(targinfoji.r,2));

            for jj = 1:size(targinfoji.r,2)
                dists = vecnorm(targinfoji.r(:,jj)-srcinfo.r(:,:));
                [dists,idst] = sort(dists);
                jjclose(jj) = idst(2);
                dclosest = min(dclosest,dists(1));
            end
            
            iside = sign(sum((targinfoji.r-srcinfo.r(:,jjclose)).*srcinfo.n(:,jjclose),1));
            iiin = iside < 0;
            iout = iside >= 0;
        end

        if any(iiin)
            targinfouse= [];
            targinfouse.r = targinfoji.r(:,iiin);
            if isfield(targinfoji, 'd')
                targinfouse.d = targinfoji.d(:,iiin);
            end
            if isfield(targinfoji, 'd2')
                targinfouse.d2 = targinfoji.d2(:,iiin);
            end
            if isfield(targinfoji, 'n')
                targinfouse.n = targinfoji.n(:,iiin);
            end  
        
            optsuse.side = 'i';

        % Helsing-Ojala (interior/exterior?)
        allmatsf = cell(size(kern.splitinfo.type));
        [allmatsf{:}] = chnk.pquadwts(r,d,n,d2,wts,i,targinfouse.r,t,w, ...
            optsuse,intp_ab,intp,kern.splitinfo.type,true);

        r_i = intp*(r(1,:,i)'+1i*r(2,:,i)'); 
        d_i = (intp*(d(1,:,i)'+1i*d(2,:,i)'));
        d2_i = (intp*(d(1,:,i)'+1i*d(2,:,i)'));
        sp = abs(d_i); tang = d_i./sp; 
        n_i = -1i*tang; 
        srcinfof = [];
        srcinfof.r  = [real(r_i)  imag(r_i)]';
        srcinfof.d  = [real(d_i)  imag(d_i)]';
        srcinfof.d2 = [real(d2_i) imag(d2_i)]';
        srcinfof.n  = [real(n_i)  imag(n_i)]';

        mat1f = zeros(opdims(1)*size(targinfouse.r,2),opdims(2)*2*k);
        funsf = kern.splitinfo.functions(srcinfof,targinfouse);
        for l = 1:length(allmatsf)
            switch kern.splitinfo.action{l}
                case 'r'
                    mat0 = real(allmatsf{l});
                case 'i'
                    mat0 = imag(allmatsf{l});
                case 'c'
                    mat0 = allmatsf{l};
            end
            mat0opdim = kron(mat0,ones(opdims(:).'));
            mat0xsplitfun = mat0opdim.*funsf{l};
            mat1f = mat1f + mat0xsplitfun;
        end
        mat1 = mat1f*kron(intp,eye(opdims(:).'));

        else
            mat1 = [];
        end
        if any(iout)
            targinfouse= [];
            targinfouse.r = targinfoji.r(:,iout);
            if isfield(targinfoji, 'd')
                targinfouse.d = targinfoji.d(:,iout);
            end
            if isfield(targinfoji, 'd2')
                targinfouse.d2 = targinfoji.d2(:,iout);
            end
            if isfield(targinfoji, 'n')
                targinfouse.n = targinfoji.n(:,iout);
            end  
        
            optsuse.side = 'e';
            % Helsing-Ojala (interior/exterior?)
            allmatsf = cell(size(kern.splitinfo.type));
            [allmatsf{:}] = chnk.pquadwts(r,d,n,d2,wts,i,targinfouse.r,t,w, ...
                optsuse,intp_ab,intp,kern.splitinfo.type,true);
    
            r_i = intp*(r(1,:,i)'+1i*r(2,:,i)'); 
            d_i = (intp*(d(1,:,i)'+1i*d(2,:,i)'));
            d2_i = (intp*(d(1,:,i)'+1i*d(2,:,i)'));
            sp = abs(d_i); tang = d_i./sp; 
            n_i = -1i*tang; 
            srcinfof = [];
            srcinfof.r  = [real(r_i)  imag(r_i)]';
            srcinfof.d  = [real(d_i)  imag(d_i)]';
            srcinfof.d2 = [real(d2_i) imag(d2_i)]';
            srcinfof.n  = [real(n_i)  imag(n_i)]';
    
            mat2f = zeros(opdims(1)*size(targinfouse.r,2),opdims(2)*2*k);
            funsf = kern.splitinfo.functions(srcinfof,targinfouse);
            for l = 1:length(allmatsf)
                switch kern.splitinfo.action{l}
                    case 'r'
                        mat0 = real(allmatsf{l});
                    case 'i'
                        mat0 = imag(allmatsf{l});
                    case 'c'
                        mat0 = allmatsf{l};
                end
                mat0opdim = kron(mat0,ones(opdims(:).'));
                mat0xsplitfun = mat0opdim.*funsf{l};
                mat2f = mat2f + mat0xsplitfun;
            end
            mat2 = mat2f*kron(intp,eye(opdims(:).'));
        else
            mat2 = [];
        end

        mat3 = zeros(size(targinfoji.r,2),k);
        mat3(iiin,:) = mat1;
        mat3(iout,:) = mat2;


        js1 = jmat:jmatend;
        js1 = repmat( (js1(:)).',opdims(1)*numel(ji),1);
        
        
        indji = (ji-1)*opdims(1);
        indji = repmat( (indji(:)).', opdims(1),1) + ( (1:opdims(1)).');
        indji = indji(:);
        
        indji = repmat(indji,1,opdims(2)*k);
        
        iend = istart+numel(mat3)-1;
        is(istart:iend) = indji(:);
        js(istart:iend) = js1(:);
        vs(istart:iend) = mat3(:);
        istart = iend+1;
    end
    mat = sparse(is,js,vs,opdims(1)*nt,opdims(2)*chnkr.npt);
    
end

if dclosest < 1e-10
    warning('Unable to estimate pquad side. Provide opts.side to ensure accuracy.')
end

end
