function obj = zeros(m, n)
%KERNEL.ZEROS   Construct the zero kernel.
%   KERNEL.ZEROS() constructs the zero kernel with operator dimensions
%   of 1 x 1.
%
%   KERNEL.ZEROS(N) constructs the zero kernel with operator dimensions
%   of N x N.
%
%   KERNEL.ZEROS(M, N) constructs the zero kernel with operator dimensions
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
        out = zeros(m*nt, n*ns);
    end

    function varargout = fmm_(eps, s, t, sigma)

        if ( isstruct(t) )
            [~, nt] = size(t.r);
        else
            [~, nt] = size(t);
        end

        if ( nargout > 0 ), varargout{1} = zeros(m*nt, 1); end
        if ( nargout > 1 ), varargout{2} = zeros(2, m*nt); end
        if ( nargout > 2 ), varargout{3} = zeros(3, m*nt); end
        if ( nargout > 3 )
            error('CHUNKIE:kernel:zeros', 'Too many output arguments for FMM.');
        end

    end

obj = kernel();
obj.name   = 'zeros';
obj.opdims = [m n];
obj.sing   = 'smooth';
obj.eval   = @eval_;
obj.fmm    = @fmm_;

end
