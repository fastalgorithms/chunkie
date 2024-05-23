function fint = chunkerintegral(chnkr,f,opts)
%CHUNKERINTEGRAL compute the integral of the function f over the chunk
% geometry. 
%
% Syntax: fint = chunkerintegral(chnkr,f,opts)
%
% Input:
%   chnkr - chunker object description of curve
%   f - either an array of values for each point on the chunker or a 
%       function handle which acts on an array of points in the
%       appropriate dimension, i.e. fvals = f(xx) is an array of values 
%       of f at each coordinate xx(1:chnkr.dim,i)
%   opts - structure for setting various parameters
%       opts.usesmooth - if = true, then just use the smooth integration
%          rule for each chunk. otherwise, adaptive integration 
%          (quadgk) is used (default 0)
%       opts.quadgkparams - if non-empty this is a cell structure
%       containing string,value pairs to be sent to quadgk (default {})
%
% output:
%   fint - integral of f over the chunker object
%
% see also QUADGK

% author: Travis Askham (askhamwhat@gmail.com)

if nargin < 3
    opts = [];
end

if ~isfield(opts,'usesmooth'); opts.usesmooth = false; end
if ~isfield(opts,'quadgkparams'); opts.quadgkparams = {}; end

k = chnkr.k;
nch = chnkr.nch;

iffun = false;
if isa(f,'function_handle')
    iffun = true;
else
    assert(isa(f,'numeric') && all(numel(f)==k*nch),...
        'f must be either a fun handle or array of size of chnkr.r');
    f = reshape(f,k,nch);
end

if opts.usesmooth
    % assume smooth weights are good enough
    wts = chnkr.wts;
    if iffun
        fvals = f(chnkr.r(:,:));
    else
        fvals = f(:);
    end
    
    fint = (wts(:)).'*fvals(:);
    
else
    % use adaptive quadrature
    
    [~,~,u] = lege.exps(k);
    
    [rc,dc] = exps(chnkr);
    
    if iffun
        fint = 0.0;
        for i = 1:nch
            rci = rc(:,:,i);
            dci = dc(:,:,i);
            fint = fint + chnk.intchunk.fhandle(f,rci,dci);
        end
    else
        fint = 0.0;
        fc = u*f;
        for i = 1:nch
            dci = dc(:,:,i);
            fci = fc(:,i);
            fint = fint + chnk.intchunk.fcoefs(fci,dci);
        end
    end
    
end
