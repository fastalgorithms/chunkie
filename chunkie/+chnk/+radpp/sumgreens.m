function g = sumgreens(coefs, varargin)
%CHNK.RADPP.SUMGREENS linear combination of Green's functions, for
% KERNEL.RADGREEN and KERNEL.RADPP.
%
% Syntax: g = chnk.radpp.sumgreens(coefs, g1, g2, ...)
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
%    K = kernel.radpp(chnk.radpp.sumgreens([1 -1], gh, gl), 4, 's');
%
% Analytically combining kernels is often more stable, see e.g.
% CHNK.HELM2D.HELMDIFFGREEN.
%
% see also KERNEL.RADGREEN, KERNEL.RADPP

if numel(coefs) ~= numel(varargin)
    error('CHNK.RADPP.SUMGREENS: need one coefficient per Green''s function.');
end
greens = varargin;
g = @(s,t) green_sum(coefs, greens, s, t);

end


function varargout = green_sum(coefs, greens, s, t)
nout = max(nargout, 1);
varargout = cell(1, nout);
for k = 1:numel(greens)
    out = cell(1, nout);
    [out{:}] = greens{k}(s, t);
    for i = 1:nout
        if k == 1
            varargout{i} = coefs(k)*out{i};
        else
            varargout{i} = varargout{i} + coefs(k)*out{i};
        end
    end
end
end
