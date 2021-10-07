function chnkr = chunkerfuncuni(fcurve,nch,cparams,pref)
%CHUNKERFUNC create a chunker object corresponding to a parameterized curve
%
% Syntax: chnkr = chunkerfunc(fcurve,cparams,pref)
%
% Input: 
%   fcurve - function handle of the form
%               [r,d,d2] = fcurve(t)
%            where r, d, d2 are size [dim,size(t)] arrays describing
%            position, first derivative, and second derivative of a curve
%            in dim dimensions parameterized by t.
%
% Optional input:
%   nch - number of chunks to use (16)
%	cparams - curve parameters structure (defaults)
%       cparams.ta = left end of t interval (0)
%       cparams.tb = right end of t interval (2*pi)
%   pref - chunkerpref object or structure (defaults)
%       pref.nchmax - maximum number of chunks (10000)
%       pref.k - number of Legendre nodes on chunks (16)
%
% Examples:
%   chnkr = chunkerfuncuni(@(t) starfish(t)); % chunk up starfish w/ standard
%                                        % options
%   pref = []; pref.k = 30; 
%   cparams = []; cparams.eps = 1e-3;
%   chnkr = chunkerfunc(@(t) starfish(t),cparams,pref); % change up options
%   
% see also CHUNKERPOLY, CHUNKERPREF, CHUNKER

% author: Travis Askham (askhamwhat@gmail.com)
%


if nargin < 2
    nch = 16;
end
if nargin < 3
    cparams = [];
end
if nargin < 4
    pref = chunkerpref();
else
    pref = chunkerpref(pref);
end



ta = 0.0; tb = 2*pi; ifclosed=true;
chsmall = Inf; nover = 0;
eps = 1.0e-6;
lvlr = 'a'; maxchunklen = Inf; lvlrfac = 2.0;

if isfield(cparams,'ta')
    ta = cparams.ta;
end	 
if isfield(cparams,'tb')
    tb = cparams.tb;
end	 
if isfield(cparams,'ifclosed')
    ifclosed = cparams.ifclosed;
end	 

chnkr = chunker(pref); % empty chunker
chnkr = chnkr.addchunk(nch);

ts = linspace(ta,tb,nch+1);

k = pref.k;
 
dim = checkcurveparam(fcurve,ta);
pref.dim = dim;
nout = 3;
out = cell(nout,1);

[xs,ws] = lege.exps(k);

%       . . . start chunking

ab = zeros(2,nch);
adjs = zeros(2,nch);
ab(1,:) = ts(1:end-1);
ab(2,:) = ts(2:end);

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


for i = 1:nch
    a=ab(1,i);
    b=ab(2,i);
    
    ts = a + (b-a)*(xs+1)/2;
    [out{:}] = fcurve(ts);
    chnkr.r(:,:,i) = reshape(out{1},dim,k);
    chnkr.d(:,:,i) = reshape(out{2},dim,k);
    chnkr.d2(:,:,i) = reshape(out{3},dim,k);
    chnkr.h(i) = (b-a)/2;
end

chnkr.adj = adjs(:,1:nch);

end
