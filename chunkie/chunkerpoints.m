function chnkr = chunkerpoints(src,opts)
%CHUNKERPOINTS create a chunker object given a set of Gauss-Legendre 
% panels discretizing a curve.
%
% Syntax: chnkr = chunkerpoints(src,opts)
%
% Input: 
%   src - If src is of class 'double', then it should be a 
%         dim x k x nch array. dim is the dimension in which the curve 
%         lives (typically 2), k is the number of Gauss-Legendre nodes 
%         per panel (typically 16), and nch is the total number of panels.
%   opts - a structure containing optional argumensts (defaults)
%       opts.ifclosed - (true) a boolean. True if the desired curve 
%           is closed, and False if not. Note that no checking is done. 
%           All this changes is the corresponding chunker adjacency info.
%
% Output:
%   chnkr - a chunker object containing the discretization of the domain
%   
% see also CHUNKERPOLY, CHUNKERPREF, CHUNKER
%

d  = [];
d2 = [];
if (strcmp(class(src),'double'))
    r = src;
elseif (strcmp(class(src),'struct'))
    if (isfield(src,'r'))
        r = src.r;
    else
        error('ERROR: missing field r in chunkerpoints');
    end
    if (isfield(src,'d') && isequal(size(r),size(src.d)))
        d = src.d;
    end
    if (isfield(src,'d2') && isequal(size(r),size(src.d)))
        d2 = src.d2;
    end
end

ifclosed = true;

if (nargin >1)
    if (isfield(opts,'ifclosed'))
        ifclosed = opts.ifclosed;
    end

end

pref = [];

dim = size(r,1);
k   = size(r,2);
nch = size(r,3); 
 
pref.dim = dim;
pref.k   = k;

[xs,~,us,vs] = lege.exps(k);
dermat = (vs*[lege.derpol(us); zeros(1,k)]).';


%       . . . start chunking

adjs = zeros(2,nch);

adjs(1,:) = (1:nch) - 1;
adjs(2,:) = (1:nch) + 1;

if ifclosed
    adjs(1,1)=nch;
    adjs(2,nch)=1;
else
    adjs(1,1)=-1;
    adjs(2,nch)=-1;
end



%       up to here, everything has been done in parameter space, [ta,tb]
%       . . . finally evaluate the k nodes on each chunk, along with 
%       derivatives and chunk lengths

chnkr = chunker(pref); % empty chunker
chnkr = chnkr.addchunk(nch);


for i = 1:nch
    
    rtmp  = squeeze(r(:,:,i));
    if (~isempty(d))
        dtmp = squeeze(d(:,:,i));
    else
        dtmp  = rtmp*dermat;
    end
    if (~isempty(d2))
        d2tmp = squeeze(d2(:,:,i));
    else
        d2tmp = dtmp*dermat;    
    end
    chnkr.rstor(:,:,i) = rtmp;
    chnkr.dstor(:,:,i) = dtmp;
    chnkr.d2stor(:,:,i) = d2tmp;
end

chnkr.adjstor(:,1:nch) = adjs(:,1:nch);

% Set normals
chnkr.nstor(:,:,1:nch) = normals(chnkr);

% Set weights
chnkr.wtsstor(:,1:nch) = weights(chnkr);

end
