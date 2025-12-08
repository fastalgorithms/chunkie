function obj = axissymlap2d(type, m)
%KERNEL.AXISSYMLAP2D   Construct the axissymmetric Laplace kernel.

if ( nargin < 1 )
    error('Missing Laplace kernel type.');
end

obj = kernel();
obj.name = 'axissymlaplace';
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.axissymlap2d.kern_modal(s, t, [0,0], 's',m);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern_modal(s, t, o, 's',m);
        obj.fmm = [];
        obj.sing = 'log';

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.axissymlap2d.kern_modal(s, t, [0,0], 'd',m);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern_modal(s, t, o, 'd',m);
        obj.fmm = [];
        obj.sing = 'smooth';

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern_modal(s, t, [0,0], 'sprime',m);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern_modal(s, t, o, 'sprime',m);
        obj.fmm = [];
        obj.sing = 'log';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.axissymlap2d.kern_modal(s, t, [0,0], 'dprime',m);
        obj.shifted_eval = @(s,t,o) chnk.axissymlap2d.kern_modal(s, t, o, 'dprime',m);
        obj.fmm = [];
        obj.sing = 'hs';

    otherwise
        error('Unknown axissym Laplace kernel type ''%s''.', type);

end

end
