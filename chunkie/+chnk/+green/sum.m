function g = sum(coefs, varargin)
%CHNK.GREEN.SUM linear combination of Green's functions or fmms, for
% KERNEL.GREEN and KERNEL.RADPP.
%
% Syntax: g = chnk.green.sum(coefs, g1, g2, ...)
%
% returns the function handle
%
%    [val, grad, hess, ...] = g(src, targ)
%
% with each output sum_k coefs(k)*(that output of gk).  
%
% Example: Helmholtz minus Laplace, fitted as one profile
%    gh = @(s,t) chnk.helm2d.green(1.1, s, t);
%    gl = @(s,t) chnk.lap2d.green(s, t);
%    K = kernel.radpp(chnk.green.sum([1 -1], gh, gl), 4, 's');
%
% Analytically combining kernels is often more stable, see e.g.
% CHNK.HELM2D.HELMDIFFGREEN.
%
% see also KERNEL.GREEN, KERNEL.RADPP

if numel(coefs) ~= numel(varargin)
    error('CHNK.GREEN.SUM: need one coefficient per Green''s function.');
end
greens = varargin;
g = @(varargin) green_sum(coefs, greens, varargin{:});

end


function varargout = green_sum(coefs, greens, varargin)
nout = max(nargout, 1);
varargout = cell(1, nout);
for k = 1:numel(greens)
    out = cell(1, nout);
    [out{:}] = greens{k}(varargin{:});
    for i = 1:nout
        varargout{i} = add(varargout{i}, coefs(k), out{i});
    end
end
end


function a = add(a, c, b)
% a + c*b, field by field for the output struct of an fmm
if ( ~isstruct(b) )
    if ( isempty(a) ), a = c*b; else, a = a + c*b; end
    return
end
for f = fieldnames(b).'
    v = c*b.(f{1});
    if ( isfield(a, f{1}) ), v = a.(f{1}) + v; end
    a.(f{1}) = v;
end
end
