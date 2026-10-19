function Kg = grad(K)
%GRAD   Target gradient of a kernel from KERNEL.GREEN, KERNEL.RADPP,
% KERNEL.LAP2D or KERNEL.HELM2D, the kernel of the gradient type, e.g.
% 'sg' for 's'.
%
% See also KERNEL.NORMALTRACE.

type = K.type;
if ( strcmp(type, 'trep') ), type = 'trans_rep'; end
if ( ~any(strcmp(type, {'s', 'd', 'c', 'trans_rep'})) )
    error('KERNEL.GRAD: Currently only supported for the types s, d, c and trans_rep.');
end
if ( strcmp(type, 'trans_rep') )
    type = 'trans_rep_grad';
else
    type = [type 'g'];
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
            Kg = kernel.green(p.green, type, args{:});
        else
            Kg = kernel.radpp(p.green, [p.rmin p.rmax], type, args{:}, ...
                'opts', p.fitopts);
        end
    case 'laplace'
        Kg = kernel.lap2d(type, coefs{:});
    case 'helmholtz'
        Kg = kernel.helm2d(type, p.zk, coefs{:});
    otherwise
        error('KERNEL.GRAD: Currently only supported for kernel.green, kernel.radpp, kernel.lap2d and kernel.helm2d.');
end

end
