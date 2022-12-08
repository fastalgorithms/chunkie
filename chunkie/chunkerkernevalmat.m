function mat = chunkerkernevalmat(chnkr,kern,targs,opts)
%CHUNKERKERNEVALMAT compute the matrix which maps density values on 
% the chunk geometry to the value of the convolution of the given
% integral kernel with the density at the specified target points
%
% Syntax: mat = chunkerkerneval(chnkr,kern,targs,opts)
%
% Input:
%   chnkr - chunker object description of curve
%   kern - integral kernel taking inputs kern(srcinfo,targinfo) 
%   targs - targ(1:2,i) gives the coords of the ith target
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

%
% output:
%   mat - (opdims(1)*nt) x opdims(2)*chnkr.npt matrix mapping values 
%         of a density on the boundary to values of convolution at target 
%         points
%
% see also CHUNKERKERNEVAL


% author: Travis Askham (askhamwhat@gmail.com)

% determine operator dimensions using first two points


srcinfo = []; targinfo = [];
srcinfo.r = chnkr.r(:,1); srcinfo.d = chnkr.d(:,1); 
srcinfo.n = chnkr.n(:,1);
srcinfo.d2 = chnkr.d2(:,1);
targinfo.r = chnkr.r(:,2); targinfo.d = chnkr.d(:,2); 
targinfo.d2 = chnkr.d2(:,2); targinfo.n = chnkr.n(:,2);

ftemp = kern(srcinfo,targinfo);
opdims = size(ftemp);

if nargin < 4
    opts = [];
end

forcesmooth = false;
forceadap = false;
forcepqud = false;
nonsmoothonly = false;
fac = 1.0;
eps = 1e-12;
if isfield(opts,'forcesmooth'); forcesmooth = opts.forcesmooth; end
if isfield(opts,'forceadap'); forceadap = opts.forceadap; end
if isfield(opts,'forcepquad'); forcepqud = opts.forcepquad; end
if isfield(opts,'nonsmoothonly'); nonsmoothonly = opts.nonsmoothonly; end
if isfield(opts,'fac'); fac = opts.fac; end
if isfield(opts,'eps'); eps = opts.eps; end

[dim,~] = size(targs);

if (dim ~= 2); warning('only dimension two tested'); end

optssmooth = []; 
optsadap = []; 
optsadap.eps = eps;

if forcesmooth
    mat = chunkerkernevalmat_smooth(chnkr,kern,opdims,targs, ...
        [],optssmooth);
    return
end

if forceadap
    mat = chunkerkernevalmat_adap(chnkr,kern,opdims, ...
        targs,[],optsadap);
    return
end

if forcepqud
    optsflag = []; optsflag.fac = fac;
    flag = flagnear(chnkr,targs,optsflag);
    spmat = chunkerkernevalmat_ho(chnkr,kern,opdims, ...
        targs,flag,optsadap);
    mat = chunkerkernevalmat_smooth(chnkr,kern,opdims,targs, ...
        flag,opts);
    mat = mat + spmat;
    return
end

% smooth for sufficiently far, adaptive otherwise

optsflag = []; optsflag.fac = fac;
flag = flagnear(chnkr,targs,optsflag);
spmat = chunkerkernevalmat_adap(chnkr,kern,opdims, ...
        targs,flag,optsadap);

if nonsmoothonly
    mat = spmat;
    return;
end

mat = chunkerkernevalmat_smooth(chnkr,kern,opdims,targs, ...
    flag,opts);

mat = mat + spmat;


end



function mat = chunkerkernevalmat_smooth(chnkr,kern,opdims, ...
    targs,flag,opts)

if nargin < 6
    flag = [];
end
if nargin < 7
    opts = [];
end

k = chnkr.k;
nch = chnkr.nch;

targinfo = []; targinfo.r = targs;
srcinfo = []; srcinfo.r = chnkr.r(:,:); srcinfo.n = chnkr.n(:,:);
srcinfo.d = chnkr.d(:,:); srcinfo.d2 = chnkr.d2(:,:);

mat = kern(srcinfo,targinfo);
wts = weights(chnkr);
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
    targs,flag,opts)

k = chnkr.k;
nch = chnkr.nch;

if nargin < 5
    flag = [];
end
if nargin < 6
    opts = [];
end

[~,nt] = size(targs);

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
    h = chnkr.h;
    targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);    
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        mat(:,jmat:jmatend) =  chnk.adapgausswts(r,d,n,d2,h,ct,bw,i,targs, ...
                    targd,targd2,kern,opdims,t,w,opts);
                
        js1 = jmat:jmatend;
        js1 = repmat( (js1(:)).',1,opdims(1)*numel(ji));
                
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
    h = chnkr.h;
    targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
    targn = zeros(chnkr.dim,nt);
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        [ji] = find(flag(:,i));
        mat1 =  chnk.adapgausswts(r,d,n,d2,h,ct,bw,i,targs(:,ji), ...
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
    targs,flag,opts)

k = chnkr.k;
nch = chnkr.nch;

if nargin < 5
    flag = [];
end
if nargin < 6
    opts = [];
end

[~,nt] = size(targs);

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
    h = chnkr.h;

    % interpolation matrix 
    intp = lege.matrin(k,t);          % interpolation from k to 2*k
    intp_ab = lege.matrin(k,[-1;1]);  % interpolation from k to end points
    targd = zeros(chnkr.dim,nt); targd2 = zeros(chnkr.dim,nt);
    targn = zeros(chnkr.dim,nt);
    for i = 1:nch
        jmat = 1 + (i-1)*k*opdims(2);
        jmatend = i*k*opdims(2);
                        
        [ji] = find(flag(:,i));

        % Helsing-Ojala (interior/exterior?)
        mat1 = chnk.pquadwts(r,d,n,d2,h,ct,bw,i,targs(:,ji), ...
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