function obj = zeros(m,n)
%KERNEL.ZEROS   Construct the zero kernel with opdims (m,n)

if ( nargin < 1 )
    error('Missing opdims');
end

if ( nargin < 2 )
    n = m;
end

obj = kernel();
obj.name = 'zeros';
obj.opdims = [m n];
obj.sing = 'smooth';

    function out = eval_(s, t)

        [~, ns] = size(s.r);
        [~, nt] = size(t.r);
        out = zeros(obj.opdims(1)*nt, obj.opdims(2)*ns);
    end

    function varargout = fmm_(eps, s, t, sigma)

        if isa(t,'struct')
            [~,nt] = size(t.r);
        else
            [~,nt] = size(t);
        end


        if ( nargout > 0 )
            varargout{1} = zeros(obj.opdims(1)*nt, 1);
        end

        if ( nargout > 1 )
            varargout{2} = zeros(2, obj.opdims(1)*nt);
        end
        
        if ( nargout > 2 )
            varargout{3} = zeros(3, obj.opdims(1)*nt);
        end

        if ( nargout > 3 )
            error('CHUNKIE:kernel:zeros', 'Too many output arguments for FMM.');
        end
    end

obj.eval = @eval_;
obj.fmm = @fmm_;

end