function obj = nans(m, n)
%KERNEL.NAN   Construct the nan kernel.
%   KERNEL.NAN() constructs the nan kernel with operator dimensions
%   of 1 x 1.
%
%   KERNEL.NAN(N) constructs the nan kernel with operator dimensions
%   of N x N.
%
%   KERNEL.NAN(M, N) constructs the NAN kernel with operator dimensions
%   of M x N.

if ( nargin < 1 )
    m = 1;
end

if ( nargin < 2 )
    n = m;
end

    function out = eval_(s, t)
        [~, ns] = size(s.r);
        [~, nt] = size(t.r);
        out = nan(m*nt, n*ns);
    end

    function out = shifted_eval_(s, t, o)
        [~, ns] = size(s.r);
        [~, nt] = size(t.r);
        out = nan(m*nt, n*ns);
    end

    function varargout = fmm_(eps, s, t, sigma)

        if ( isstruct(t) )
            [~, nt] = size(t.r);
        else
            [~, nt] = size(t);
        end

        if ( nargout > 0 ), varargout{1} = nan(m*nt, 1); end
        if ( nargout > 1 ), varargout{2} = nan(2, m*nt); end
        if ( nargout > 2 ), varargout{3} = nan(3, m*nt); end
        if ( nargout > 3 )
            error('CHUNKIE:kernel:nan', 'Too many output arguments for FMM.');
        end

    end

obj = kernel();
obj.name   = 'nans';
obj.opdims = [m n];
obj.sing   = 'smooth';
obj.eval   = @eval_;
obj.shifted_eval   = @shifted_eval_;
obj.fmm    = @fmm_;
obj.isnan = true;

end
