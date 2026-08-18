function [x,w,u,v] = exps(k,opts)
%LEGE.EXPS get nodes, weights, and matrices for switching between values
% and coefficients for Legendre nodes of the specified order
%
% IN:
%   k - the desired order (number of points)
%   opts - options structure. 
%       opts.stab (false): if true, force use of slow but more stable
%       method for weights and nodes. default false uses linear scaling
%       routine which can lose some precision in the weights
% OUT: 
%   x - k Legendre nodes on [-1,1]
%   w - the corresponding integration weights
%   u - matrix which maps function values at k Legendre nodes to the 
%   corresponding coefficients
%   v - matrix which evaluates the Legendre polynomial with the given 
%   coefficients at the k Legendre nodes.
% 
stab = false;
if nargin > 1
    if isfield(opts,'stab')
        stab = opts.stab;
    end
end

% Cache outputs keyed by (k,stab,nargout); kernel-split / RCIP loops
% hit this thousands of times per build.
persistent cache
if isempty(cache); cache = struct(); end
key = sprintf('k%d_s%d_n%d', k, stab, nargout);
if isfield(cache, key)
    e = cache.(key);
    x = e.x;
    if nargout > 1; w = e.w; end
    if nargout > 2; u = e.u; end
    if nargout > 3; v = e.v; end
    return
end

if (or(stab,k<=200))
    if nargout > 1
        [x,w] = lege.rts_stab(k);
    else
        x = lege.rts_stab(k);
    end
else
    if nargout > 1
        [x,w] = lege.rts(k);
    else
        x = lege.rts(k);
    end
end
    
if nargout > 2
    v = (lege.pols(x(:),k-1)).';
    d = (2.0*(1:k) - 1)/2.0;
    u = ((v).*(w(:)*d)).';
end

e = struct('x', x);
if nargout > 1; e.w = w; end
if nargout > 2; e.u = u; end
if nargout > 3; e.v = v; end
cache.(key) = e;
