function obj = lap2d(type, coefs)
%KERNEL.LAP2D   Construct the Laplace kernel.
%   KERNEL.LAP2D('s') or KERNEL.LAP2D('single') constructs the single-layer
%   Laplace kernel.
%
%   KERNEL.LAP2D('d') or KERNEL.LAP2D('double') constructs the double-layer
%   Laplace kernel.
%
%   KERNEL.LAP2D('sp') or KERNEL.LAP2D('sprime') constructs the derivative
%   of the single-layer Laplace kernel.
%
%   KERNEL.LAP2D('c', coefs) or KERNEL.LAP2D('combined', coefs) constructs
%   the combined-layer Laplace kernel with parameter ETA, i.e.,
%   KERNEL.LAP2D('d') + ETA*KERNEL.LAP2D('s').
%
% See also CHNK.LAP2D.KERN.

if ( nargin < 1 )
    error('Missing Laplace kernel type.');
end

obj = kernel();
obj.name = 'laplace';
obj.opdims = [1 1];

switch lower(type)

    case {'s', 'single'}
        obj.type = 's';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 's');
        obj.fmm  = @(eps,s,t,sigma) chnk.lap2d.fmm(eps, s, t, 's', sigma);
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[1 0 0 0]};
        obj.splitinfo.action = {'r'};
        obj.splitinfo.function = {@(s,t) ones(size(t.r,2),size(s.r,2))};

    case {'d', 'double'}
        obj.type = 'd';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'd');
        obj.fmm  = @(eps,s,t,sigma) chnk.lap2d.fmm(eps, s, t, 'd', sigma);
        obj.sing = 'smooth';
        obj.splitinfo = [];
        obj.splitinfo.type = {[0 0 -1 0]};
        obj.splitinfo.action = {'r'};
        obj.splitinfo.function = {@(s,t) ones(size(t.r,2),size(s.r,2))};

    case {'sp', 'sprime'}
        obj.type = 'sp';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'sprime');
        obj.fmm  = @(eps,s,t,sigma) chnk.lap.fmm(eps, s, t, 'sprime', sigma);
        obj.sing = 'smooth';

    case {'st', 'stau'}
        obj.type = 'st';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'stau');
        obj.fmm  = @(eps,s,t,sigma) chnk.lap.fmm(eps, s, t, 'stau', sigma);
        obj.sing = 'pv';

    case {'dp', 'dprime'}
        obj.type = 'dp';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'dprime');
        obj.fmm  = @(eps,s,t,sigma) chnk.lap.fmm(eps, s, t, 'dprime', sigma);
        obj.sing = 'hs';

    case {'sg', 'sgrad'}
        obj.type = 'sg';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'sgrad');
        obj.fmm = @(eps,s,t,sigma) chnk.lap2d.fmm(eps, s, t, 'sgrad', sigma);
        obj.sing = 'pv';
        obj.opdims = [2,1];

    case {'dg', 'dgrad'}
        obj.type = 'dg';
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'dgrad');
        obj.fmm = @(eps,s,t,sigma) chnk.lap2d.fmm(eps, s, t, 'dgrad', sigma);
        obj.sing = 'hs';
        obj.opdims = [2,1];

    case {'c', 'combined'}
        if ( nargin < 2 )
            warning('Missing combined layer parameter coefs. Defaulting to 1.');
            coefs = ones(2,1);
        end
        obj.type = 'c';
        obj.params.coefs = coefs;
        obj.eval = @(s,t) chnk.lap2d.kern(s, t, 'c', coefs);
        obj.fmm  = @(eps,s,t,sigma) chnk.lap2d.fmm(eps, s, t, 'c', sigma, coefs);
        obj.sing = 'log';
        obj.splitinfo = [];
        obj.splitinfo.type = {[1 0 0 0],[0 0 -1 0]};
        obj.splitinfo.action = {'r','r'};
        obj.splitinfo.function = cell(2,1);
        obj.splitinfo.function{1} = @(s,t) coefs(2)*ones(size(t.r,2),size(s.r,2));
        obj.splitinfo.function{2} = @(s,t) coefs(2)*ones(size(t.r,2),size(s.r,2));

    otherwise
        error('Unknown Laplace kernel type ''%s''.', type);

end

icheck = exist(['fmm2d.' mexext], 'file');
if icheck ~=3
    obj.fmm = [];
end

end
