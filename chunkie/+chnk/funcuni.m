function chnkr = funcuni(fcurve,cparams,pref)
%CHNK.FUNCUNI uniform discretization of parameterized curve, mostly for
%   debugging and convergence tests.
%
% Syntax: chnkr = chnk.funcuni(fcurve,cparams,pref)
%
% Input: 
%   fcurve - function handle of the form
%               [r,d,d2] = fcurve(t)
%            where r, d, d2 are size [dim,size(t)] arrays describing
%            position, first derivative, and second derivative of a curve
%            in dim dimensions parameterized by t.
%
% Optional input:
%	cparams - curve parameters structure (defaults)
%   	cparams.ta = left end of t interval (0)
%       cparams.tb = right end of t interval (2*pi)
%       cparams.npan = number of panels
%   pref - chunkerpref object or structure (defaults)
%       pref.nchmax - maximum number of chunks (10000)
%       pref.k - number of Legendre nodes on chunks (16)
%
% Examples:
%   chnkr = chnk.funcuni(@(t) starfish(t)); % chunk up starfish w/ standard
%                                        % options
%   pref = []; pref.k = 30; 
%   cparams = []; cparams.npan = 30;
%   chnkr = chnk.funcuni(@(t) starfish(t),cparams,pref); 
%                                           % change up options
%   
% see also CHUNKERFUNC, CHUNKERPOLY, CHUNKERPREF, CHUNKER

% author: Travis Askham (askhamwhat@gmail.com)
%

if nargin < 2
    cparams = [];
end
if nargin < 3
    pref = chunkerpref();
else
    pref = chunkerpref(pref);
end


ta = 0.0; tb = 2*pi; 
npan = 16;

if isfield(cparams,'ta')
    ta = cparams.ta;
end	 
if isfield(cparams,'tb')
    tb = cparams.tb;
end	 
if isfield(cparams,'npan')
    npan = cparams.npan;
end

k = pref.k;
nchmax = pref.nchmax; 
 
dim = checkcurveparam(fcurve,ta);
pref.dim = dim;
nout = 3;
out = cell(nout,1);


ifprocess = zeros(nchmax,1);

%       construct legendre nodes and weights, k and 2k of them, as well
%       as the interpolation/coefficients matrices

xs = lege.exps(k);

%       . . . start chunking

nch = npan;

ab = zeros(2,nch);
adjs = zeros(2,nch);
ab(1,1)=ta;
ab(2,1)=tb;

h = (tb-ta)/nch;

for i = 1:nch
    ab(1,i) = ta+(i-1)*h;
    ab(2,i) = ta+i*h;    
    adjs(1,i) = i-1;
    adjs(2,i) = i+1;
end
    
adjs(1,1) = nch;
adjs(2,nch) = 1;


%       up to here, everything has been done in parameter space, [ta,tb]
%       . . . finally evaluate the k nodes on each chunk, along with 
%       derivatives and chunk lengths

chnkr = chunker(pref); % empty chunker
chnkr = chnkr.addchunk(nch);

for i = 1:nch
    a=ab(1,i);
    b=ab(2,i);
    h = (b-a)/2;
    
    ts = a + (b-a)*(xs+1)/2;
    [out{:}] = fcurve(ts);
    chnkr.r(:,:,i) = reshape(out{1},dim,k);
    chnkr.d(:,:,i) = reshape(out{2},dim,k)*h;
    chnkr.d2(:,:,i) = reshape(out{3},dim,k)*h*h;
end

chnkr.adj = adjs(:,1:nch);
chnkr.wtsstor(:,1:nch) = weights(chnkr);
chnkr.nstor(:,:,1:nch) = normals(chnkr);

end


