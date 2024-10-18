    function mat = chunkerkernevalmat(chnkr,kern,targobj,opts)
%CHUNKERKERNEVALMAT compute the matrix which maps density values on 
% the chunk geometry to the value of the convolution of the given
% integral kernel with the density at the specified target points
%
% Syntax: mat = chunkerkernevalmat(chnkr,kern,targs,opts)
%
% Input:
%   chnkr - chunker object describing boundary, currently
%              only supports chunkers, and not chunkgraphs
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
%                         (e.g. self and neighbor interactoins) and return
%                         in a sparse array.
%           opts.corrections = boolean (false), if true, only compute the
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

% convert kernel to kernel object, put in singularity info 
% opts.sing provides a default value for singularities if not 
% defined for kernels

if isa(kern,'function_handle')
    kern2 = kernel(kern);
    kern = kern2;
elseif isa(kern,'cell')
    sz = size(kern);
    kern2(sz(1),sz(2)) = kernel();
    for j = 1:sz(2)
        for i = 1:sz(1)
            if isa(kern{i,j},'function_handle')
                kern2(i,j) = kernel(kern{i,j});
            elseif isa(kern{i,j},'kernel')
                kern2(i,j) = kern{i,j};
            else
                msg = "Second input is not a kernel object, function handle, " ...
                    + "or cell array";
                error(msg);
            end
        end
    end
    kern = kern2;
    
elseif ~isa(kern,'kernel')
    msg = "Second input is not a kernel object, function handle, " ...
                + "or cell array";
    error(msg);
end
    
% determine operator dimensions using first two points


srcinfo = []; targinfo = [];
srcinfo.r = chnkr.r(:,1); srcinfo.d = chnkr.d(:,1); 
srcinfo.n = chnkr.n(:,1);
srcinfo.d2 = chnkr.d2(:,1);
targinfo.r = chnkr.r(:,2); targinfo.d = chnkr.d(:,2); 
targinfo.d2 = chnkr.d2(:,2); targinfo.n = chnkr.n(:,2);

ftmp = kern.eval;

ftemp = ftmp(srcinfo,targinfo);
opdims = size(ftemp);

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
elseif isa(targobj, "chunkgraph")
    targinfo.r = targobj.r(:,:);
    targinfo.d = targobj.d(:,:);
    targinfo.d2 = targobj.d2(:,:);
    targinfo.n = targobj.n(:,:);
else
    targinfo.r = targobj;
end

[dim,~] = size(targinfo.r);


if (dim ~= 2); warning('only dimension two tested'); end

optssmooth = []; 
optsadap = []; 
optsadap.eps = eps;



if forcesmooth
    mat = chunkerkernevalmat_smooth(chnkr,ftmp,opdims,targinfo, ...
        [],optssmooth);
    return
end


if corrections
    mat = chunkerkernevalmat_adap(chnkr,ftmp,opdims, ...
        targinfo,[],optsadap);
    mat = mat-chunkerkernevalmat_smooth(chnkr,ftmp,opdims,targinfo, ...
        [],opts);
    mat = sparse(mat);
    return
end

if forceadap
    mat = chunkerkernevalmat_adap(chnkr,ftmp,opdims, ...
        targinfo,[],optsadap);
    return
end



if forcepquad
    optsflag = []; optsflag.fac = fac;
    flag = flagnear(chnkr,targinfo.r,optsflag);
    spmat = chunkerkernevalmat_ho(chnkr,ftmp,opdims, ...
        targinfo,flag,optsadap);
    mat = chunkerkernevalmat_smooth(chnkr,ftmp,opdims,targinfo, ...
        flag,opts);
    mat = mat + spmat;
    return
end

% smooth for sufficiently far, adaptive otherwise

% TODO: change to chunkerkerneval system, need routine to generate
% upsampling matrix.

optsflag = []; optsflag.fac = fac;
flag = flagnear(chnkr,targinfo.r,optsflag);
spmat = chunkerkernevalmat_adap(chnkr,ftmp,opdims, ...
        targinfo,flag,optsadap);

if nonsmoothonly
    mat = spmat;
    return;
end

mat = chunkerkernevalmat_smooth(chnkr,ftmp,opdims,targinfo, ...
    flag,opts);

mat = mat + spmat;


end



function mat = chunkerkernevalmat_smooth(chnkr,kern,opdims, ...
    targinfo,flag,opts)

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

mat = kern(srcinfo,targinfo);
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

function mat = chunkerkernevalmat_adap(chnkr,kern,opdims, ...
    targinfo,flag,opts)

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
    
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        mat(:,jmat:jmatend) =  chnk.adapgausswts(r,d,n,d2,ct,bw,i,targs, ...
                    targd,targn,targd2,kern,opdims,t,w,opts);
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
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        [ji] = find(flag(:,i));
        mat1 =  chnk.adapgausswts(r,d,n,d2,ct,bw,i,targs(:,ji), ...
                    targd(:,ji),targn(:,ji),targd2(:,ji),kern,opdims,t,w,opts);
                
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

function mat = chunkerkernevalmat_ho(chnkr,kern,opdims, ...
    targinfo,flag,opts)

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

    % interpolation matrix 
    intp = lege.matrin(k,t);          % interpolation from k to 2*k
    intp_ab = lege.matrin(k,[-1;1]);  % interpolation from k to end points
    
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        [ji] = find(flag(:,i));

        % Helsing-Ojala (interior/exterior?)
        mat1 = chnk.pquadwts(r,d,n,d2,ct,bw,i,targs(:,ji), ...
                    targd(:,ji),targn(:,ji),targd2(:,ji),kern,opdims,t,w,opts,intp_ab,intp); % depends on kern, different mat1?
                
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
