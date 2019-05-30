function [x,w,u,v] = exps(k,opts)

stab = false;
if nargin > 1
    if isfield(opts,'stab')
        stab = opts.stab;
    end
end

if (stab)
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
