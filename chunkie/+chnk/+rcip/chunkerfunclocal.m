function chnkr = chunkerfunclocal(fcurve,ts,pref,xs,ws)
%CHUNKERFUNC create a chunker object for constructing the system matrix
% used in the forward recursion for computing the preconditioner R 
% in the RCIP method
% 
%
% Input: 
%   fcurve - function handle of the form
%               [r,d,d2] = fcurve(t)
%            where r, d, d2 are size [dim,size(t)] arrays describing
%            position, first derivative, and second derivative of a curve
%            in dim dimensions parameterized by t.
%
% Optional input:
%   ts - endpoints in the parameter space
%
%   pref - chunkerpref object or structure (defaults)
%       pref.nchmax - maximum number of chunks (10000)
%       pref.k - number of Legendre nodes on chunks (16)
%   xs - Gauss-Legendre nodes of order k
%   ws - Gauss-Legendre weights of order k

if nargin < 3 || isempty(pref)
    pref = chunkerpref();
else
    pref = chunkerpref(pref);
end

if nargin < 5
    [xs,ws] = lege.exps(pref.k);
end
    

chnkr = chunker(pref,xs,ws); % empty chunker

nch = length(ts)-1; % number of chunks

chnkr = chnkr.addchunk(nch);

k = pref.k;
 
dim = checkcurveparam(fcurve,ts(1));
pref.dim = dim;
nout = 3;
out = cell(nout,1);

%[xs,~] = lege.exps(k);

%       . . . start chunking

ab = zeros(2,nch);
adjs = zeros(2,nch);
ab(1,:) = ts(1:end-1);
ab(2,:) = ts(2:end);

adjs(1,:) = (1:nch) - 1;
adjs(2,:) = (1:nch) + 1;


adjs(1,1)=-1;
adjs(2,nch)=-1;

%       up to here, everything has been done in parameter space, [ta,tb]
%       . . . finally evaluate the k nodes on each chunk, along with 
%       derivatives and chunk lengths


for i = 1:nch
    a=ab(1,i);
    b=ab(2,i);
    
    ts = a + (b-a)*(xs+1)/2;
    [out{:}] = fcurve(ts);
    h = (b-a)/2;
    chnkr.rstor(:,:,i) = reshape(out{1},dim,k);
    chnkr.dstor(:,:,i) = reshape(out{2},dim,k)*h;
    chnkr.d2stor(:,:,i) = reshape(out{3},dim,k)*h*h;

end

chnkr.adjstor(:,1:nch) = adjs(:,1:nch);
% added by Shidong Jiang
chnkr.nstor(:,:,1:nch) = normals(chnkr);

% add weights
chnkr.wtsstor(:,1:nch) = weights(chnkr);
end
