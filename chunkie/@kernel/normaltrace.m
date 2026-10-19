function Kn = normaltrace(K)
%NORMALTRACE   Target normal derivative of a kernel from KERNEL.GREEN,
% KERNEL.RADPP, KERNEL.LAP2D, KERNEL.HELM2D or KERNEL.HELM2DDIFF, the
% kernel of the prime type, e.g. 'sp' for 's'.
%
% See also KERNEL.GRAD, KERNEL.NTRACE.

type = K.type;
if ( strcmp(type, 'trep') ), type = 'trans_rep'; end
if ( ~any(strcmp(type, {'s', 'd', 'c', 'trans_rep'})) )
    error('KERNEL.NORMALTRACE: Currently only supported for the types s, d, c and trans_rep.');
end
if ( strcmp(type, 'trans_rep') )
    type = 'trans_rep_prime';
else
    type = [type 'p'];
end

p = K.params;
coefs = {};
if ( isfield(p, 'coefs') ), coefs = {p.coefs}; end

switch K.name
    case {'green', 'radpp'}
        % remove a preset 'sing' from p.args
        args = p.args;
        i = find(cellfun(@(a) (ischar(a) || isstring(a)) && strcmpi(a, 'sing'), args));
        args([i, i+1]) = [];
        if ( strcmp(K.name, 'green') )
            Kn = kernel.green(p.green, type, args{:});
        else
            Kn = kernel.radpp(p.green, [p.rmin p.rmax], type, args{:}, ...
                'opts', p.fitopts);
        end
    case 'laplace'
        Kn = kernel.lap2d(type, coefs{:});
    case 'helmholtz'
        if ( isfield(p, 'zks') )
            Kn = kernel.helm2ddiff(type, p.zks, coefs{:});
        else
            Kn = kernel.helm2d(type, p.zk, coefs{:});
        end
    otherwise
        error('KERNEL.NORMALTRACE: Currently only supported for kernel.green, kernel.radpp, kernel.lap2d, kernel.helm2d and kernel.helm2ddiff.');
end

end
