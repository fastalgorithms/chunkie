function [mat,maxrecs,numints,iers] = adapgausswts(r,d,n,d2,data,ct,bw,j,...
    rt,dt,nt,d2t,datat,kern,opdims,t,w,opts)
%CHNK.ADAPGAUSSWTS adaptive integration for interaction of kernel on chunk 
% at targets
%
% WARNING: this routine is not designed to be user-callable and assumes 
%   a lot of precomputed values as input
% WARNING: this routine currently assumes that the kernel function is
%   translation invariant to recenter (improves stability). 
%
% Syntax: [mat,maxrecs,numints,iers] = adapgausswts(r,d,d2,ct,bw,j, ...
%   rt,dt,d2t,kern,opdims,t,w,opts)
%
% Input:
%   r - chnkr nodes
%   d - chnkr derivatives at nodes
%   n - chnkr normals at nodes
%   d2 - chnkr 2nd derivatives at nodes
%   ct - Legendre nodes at order of chunker
%   bw - barycentric interpolation weights for Legendre nodes at order of
%   chunker
%   j - chunk of interest
%   rt,dt,nt,d2t - position, derivative, normal,second derivative of select 
%               target points. if any are not used by kernel (or not well
%               defined, e.g. when not on curve), a dummy array
%               of the appropriate size should be supplied
%   kern - kernel function of form kern(srcinfo,targinfo)
%   opdims - dimensions of kernel
%   t - (Legendre) integration nodes for adaptive integration
%   w - integration nodes for adaptive integrator (t and w not necessarily 
%       same order as chunker order)
%
% Output
%   mat - integration matrix
%   maxrecs - maximum recursion depth
%   numints - number of integrals computed
%   iers - error flags per target
%

eps = 1e-12;
nnmax=100000;
maxdepth=52; % since double-precision only has 52-bits of precision
transinv = true;

if nargin < 18
    opts = [];
end

if isfield(opts,'eps')
    eps = opts.eps;
end
if isfield(opts,'maxints')
    nnmax = opts.maxints;
end
if isfield(opts,'maxdepth')
    maxdepth = opts.maxdepth;
end
if isfield(opts,'transinv')
  transinv = opts.transinv;
end

[~,ntarg] = size(rt);
[~,k,~] = size(r);
k2 = length(t);

rs_ = r(:,:,j);
ds = d(:,:,j);
ns = n(:,:,j);
d2s = d2(:,:,j);
datas = [];
if ~isempty(data)
    datas = data(:,:,j);
end

stack = zeros(2,maxdepth);
vals = zeros(opdims(1)*opdims(2)*k,maxdepth);

mat = zeros(opdims(1)*ntarg,opdims(2)*k);

numints = zeros(ntarg,1); iers = zeros(ntarg,1); maxrecs = zeros(ntarg,1);

% outer loop --- over targets

for ii = 1:ntarg
    
    rt1 = rt(:,ii);
    dt1 = dt(:,ii);
    nt1 = nt(:,ii);
    d2t1 = d2t(:,ii);
    datat1 = [];
    if ~isempty(datat)
        datat1 = datat(:,ii);
    end

    % shift everything so that the target is at the origin
    if transinv
        rs = bsxfun(@minus, rs_, rt1);
        rt1 = rt1*0;
    else
        rs = rs_;
    end

    % start the recursion

    stack(1,1)=-1;
    stack(2,1)=1;
    vals(:,1) = oneintp(-1,1,rs,ds,ns,d2s,datas,ct,bw,rt1,dt1,nt1,d2t1,datat1,kern,opdims,t,w);

    % recursively integrate the thing

    jj=1;
    mat1=0;
    iers(ii)=0;
    maxrecs(ii)=0;
    toomanyints = true;
    for i=1:nnmax
        numints(ii)=i;
        if(jj > maxrecs(ii)); maxrecs(ii)=jj; end

%       subdivide the current subinterval

        a = stack(1,jj); b = stack(2,jj);
        c=(a+b)/2;
        v2 = oneintp(a,c,rs,ds,ns,d2s,datas,ct,bw,rt1,dt1,nt1,d2t1,datat1,kern,opdims,t,w);
        v3 = oneintp(c,b,rs,ds,ns,d2s,datas,ct,bw,rt1,dt1,nt1,d2t1,datat1,kern,opdims,t,w);
    
        dd= max(abs(v2+v3-vals(:,jj)));
        if(dd <= eps) 


    %       if the function on this subinterval has been 
    %       integrated with sufficient accuracy - add the 
    %       value to that of the global integral and move up
    %       in the stack

    %
            mat1=mat1+v2+v3;
            jj=jj-1;
    %
    %        if the whole thing has been integrated - return
    %
            if(jj == 0); toomanyints = false; break; end

        else

    %       if the function on this subinterval has not been 
    %       integrated with sufficient accuracy - move 
    %       down the stack

            stack(1,jj+1)=stack(1,jj);
            stack(2,jj+1)=(stack(1,jj)+stack(2,jj))/2;
            vals(:,jj+1)=v2;
    %
            stack(1,jj)=(stack(1,jj)+stack(2,jj))/2;
            vals(:,jj)=v3;
        %
            jj=jj+1;
        %     
        %       if the depth of the recursion has become excessive - bomb
        %
            if(jj > maxdepth) 
                toomanyints = false;
                iers(ii) = 8;
                break;
            end 
        end
    end

    if toomanyints; iers(ii) = 16; end
    
    istart = (ii-1)*opdims(1)+1; iend = ii*opdims(1);
    mat(istart:iend,:) = reshape(mat1,opdims(1),opdims(2)*k);
    
end

end

function val = oneintp(a,b,rs,ds,ns,d2s,datas,ct,bw,rt,dt,nt,d2t,datat,kern,opdims,t,w)
%       integrate the kernel multiplied by each Lagrange interpolant
%   on the interval [a,b] at a single target

u=(b-a)/2;
v=(b+a)/2;
tt = u*t+v;

ntt = length(t);

interpmat = bsxfun(@rdivide,bw(:),bsxfun(@minus,ct(:),(tt(:)).'));
interpmatsum = sum(interpmat,1);
%lagrange interp values at nodes on [a,b]
interpmat = bsxfun(@rdivide,interpmat,interpmatsum); 

rint = rs*interpmat;
dint = ds*interpmat;
nint = ns*interpmat;
d2int = d2s*interpmat;
dataint = [];
if ~isempty(datas)
    dataint = datas*interpmat;
end
dintlen = sqrt(sum(dint.^2,1));
%tauint = bsxfun(@rdivide,dint,dintlen);
srcinfo = []; srcinfo.r = rint; srcinfo.d = dint; 
srcinfo.d2 = d2int; srcinfo.n = nint; srcinfo.data = dataint;
targinfo = []; targinfo.r = rt; targinfo.d = dt; 
targinfo.d2 = d2t; targinfo.n = nt; targinfo.data = datat;
mat_tt = kern(srcinfo,targinfo);

dsdt = u*( (w(:).' ).*dintlen);
dsdt = repmat(dsdt,opdims(2),1); dsdt = dsdt(:);

mat_tt = bsxfun(@times,mat_tt,dsdt.');

mat_tt = reshape(mat_tt,opdims(1)*opdims(2),ntt);
val = mat_tt*(interpmat.');
val = val(:);

end

