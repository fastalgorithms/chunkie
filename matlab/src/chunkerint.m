function fint = chunkerint(chunker,f,opts)
%CHUNKERINT compute the integral of the function f over the chunk
% geometry. 
%
% input:
%   chunker - chunks description of curve
%   f - either an array of values for each point on the chunker or a 
%       function handle which acts on an array of 2d points, 
%       i.e. fvals = f(xx) is an array of values for each xx(1:2,i) pair
%   opts - structure for setting various parameters
%       opts.usesmooth - if = 1, then just use the smooth integration
%          rule for each chunk. otherwise, adaptive integration 
%          (quadgk) is used (default 0)
%       opts.quadgkparams - if non-empty this is a cell structure
%       containing string,value pairs to be sent to quadgk (default {})
%
% output:
%   fint - 
%
% see also QUADGK

if nargin < 3
    opts = [];
end

if ~isfield(opts,'usesmooth'); opts.usesmooth = false; end
if ~isfield(opts,'quadgkparams'); opts.quadgkparams = {}; end

k = chunker.k;
nch = chunker.nch;

iffun = false;
if isa(f,'function_handle')
    iffun = true;
else
    assert(isa(f,'numeric') && all(numel(f)==k*nch),...
        'f must be either a fun handle or array of size of chunker.chunks');
    f = reshape(f,k,nch);
end

if opts.usesmooth
    % assume smooth weights are good enough
    whts = chunkwhts(chunker);
    if iffun
        fint = 0.0;
        fvals = f(reshape(chunker.chunks,2,chunker.k*chunker.nch));
    else
        fvals = f(:);
    end
    
    fint = (whts(:)).'*fvals(:);
    
else
    % use adaptive quadrature
    [xc,yc,xpc,ypc] = chunkerexps(chunker);
    [x,w,u,v] = legeexps(k);
    
    if iffun
        fint = 0.0;
        for i = 1:nch
            xci = xc(:,i); yci = yc(:,i); 
            xpci = xpc(:,i); ypci = ypc(:,i);
            fint = fint + ...
                chunkerintchunk_fhandle(f,xci,yci,xpci,ypci)*chunker.hs(i);
        end
    else
        fint = 0.0;
        fc = u*f;
        for i = 1:nch
            xci = xc(:,i); yci = yc(:,i); 
            xpci = xpc(:,i); ypci = ypc(:,i);
            fci = fc(:,i);
            fint = fint + chunkerintchunk_fcoefs(fci,xpci,ypci)*chunker.hs(i);
        end
    end
    
end