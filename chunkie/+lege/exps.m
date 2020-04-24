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
