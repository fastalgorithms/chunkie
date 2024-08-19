function varargout = kern(mu, srcinfo, targinfo, type, varargin)
%CHNK.STOK2D.KERN   Evaluate the Stokes kernel.

rx = targinfo.r(1,:).' - srcinfo.r(1,:);
ry = targinfo.r(2,:).' - srcinfo.r(2,:);
r2 = rx.^2 + ry.^2;
[Nt, Ns] = size(rx);

switch lower(type)

    case {'svel', 's'}

        r = sqrt(r2);
        log1r = log(1./r);
        Kxx = 1/(4*pi*mu) * (rx.^2 ./ r2 + log1r);
        Kyy = 1/(4*pi*mu) * (ry.^2 ./ r2 + log1r);
        Kxy = 1/(4*pi*mu) * rx.*ry ./ r2;

        if ( nargout == 1 )
            % Interleave
            K = zeros(2*Nt, 2*Ns);
            K(1:2:end, 1:2:end) = Kxx;
            K(1:2:end, 2:2:end) = Kxy;
            K(2:2:end, 1:2:end) = Kxy;
            K(2:2:end, 2:2:end) = Kyy;
            varargout = {K};
        else
            varargout = {Kxx, Kxy, Kyy};
        end

    case 'spres'

        Kx = 1/(2*pi) * rx ./ r2;
        Ky = 1/(2*pi) * ry ./ r2;

        if ( nargout == 1 )
            % Interleave
            K = zeros(Nt, 2*Ns);
            K(:,1:2:end) = Kx;
            K(:,2:2:end) = Ky;
            varargout = {K};
        else
            varargout = {Kx, Ky};
        end

    case 'strac'

        r4 = r2.^2;
        nx = targinfo.n(1,:).';
        ny = targinfo.n(2,:).';
        rn = rx.*nx + ry.*ny;
        Kxx = -1/pi * rx.^2  .* rn ./ r4;
        Kyy = -1/pi * ry.^2  .* rn ./ r4;
        Kxy = -1/pi * rx.*ry .* rn ./ r4;

        if ( nargout == 1 )
            % Interleave
            K = zeros(2*Nt, 2*Ns);
            K(1:2:end, 1:2:end) = Kxx;
            K(1:2:end, 2:2:end) = Kxy;
            K(2:2:end, 1:2:end) = Kxy;
            K(2:2:end, 2:2:end) = Kyy;
            varargout = {K};
        else
            varargout = {Kxx, Kxy, Kyy};
        end

    case {'dvel', 'd'}

        r4 = r2.^2;
        nx = srcinfo.n(1,:);
        ny = srcinfo.n(2,:);
        rn = rx.*nx + ry.*ny;
        Kxx = 1/pi * rx.^2  .* rn ./ r4;
        Kyy = 1/pi * ry.^2  .* rn ./ r4;
        Kxy = 1/pi * rx.*ry .* rn ./ r4;

        if ( nargout == 1 )
            % Interleave
            K = zeros(2*Nt, 2*Ns);
            K(1:2:end, 1:2:end) = Kxx;
            K(2:2:end, 1:2:end) = Kxy;
            K(1:2:end, 2:2:end) = Kxy;
            K(2:2:end, 2:2:end) = Kyy;
            varargout = {K};
        else
            varargout = {Kxx, Kxy, Kyy};
        end

    case 'dpres'

        r4 = r2.^2;
        nx = srcinfo.n(1,:);
        ny = srcinfo.n(2,:);
        rn = rx.*nx + ry.*ny;
        Kx = mu/pi * (-nx ./ r2 + 2*rx.*rn ./ r4);
        Ky = mu/pi * (-ny ./ r2 + 2*ry.*rn ./ r4);

        if ( nargout == 1 )
            % Interleave
            K = zeros(Nt, 2*Ns);
            K(:,1:2:end) = Kx;
            K(:,2:2:end) = Ky;
            varargout = {K};
        else
            varargout = {Kx, Ky};
        end

    case 'dtrac'

        r4 = r2.^2;
        r6 = r2.^3;
        nx_s = srcinfo.n(1,:);
        ny_s = srcinfo.n(2,:);
        nx_t = targinfo.n(1,:).';
        ny_t = targinfo.n(2,:).';
        rn_s = rx.*nx_s + ry.*ny_s;
        rn_t = rx.*nx_t + ry.*ny_t;
        nn = nx_t.*nx_s + ny_t.*ny_s;

        Kxx = mu/pi * ( -8*rx.^2.*rn_t.*rn_s ./ r6 + ...
            (rx.*nx_t.*rn_s + rx.^2.*nn + rn_t.*rn_s + nx_s.*rx.*rn_t) ./ r4 + ...
            nx_t.*ny_t ./ r2 );
        Kyy = mu/pi * ( -8*ry.^2.*rn_t.*rn_s ./ r6 + ...
            (ry.*ny_t.*rn_s + ry.^2.*nn + rn_t.*rn_s + ny_s.*ry.*rn_t) ./ r4 + ...
            ny_t.*ny_s ./ r2 );
        Kxy = mu/pi * ( -8*rx.*ry.*rn_t.*rn_s ./ r6 + ...
            (rx.*ny_t.*rn_s + rx.*ry.*nn + nx_s.*ry.*rn_t) ./ r4 + ...
            nx_t.*ny_s ./ r2 );
        Kyx = mu/pi * ( -8*ry.*rx.*rn_t.*rn_s ./ r6  + ...
            (ry.*nx_t.*rn_s + ry.*rx.*nn + ny_s.*rx.*rn_t) ./ r4 + ...
            ny_t.*nx_s ./ r2 );

        if ( nargout == 1 )
            % Interleave
            K = zeros(2*Nt, 2*Ns);
            K(1:2:end, 1:2:end) = Kxx;
            K(1:2:end, 2:2:end) = Kxy;
            K(2:2:end, 1:2:end) = Kyx;
            K(2:2:end, 2:2:end) = Kyy;
            varargout = {K};
        else
            varargout = {Kxx, Kxy, Kyy};
        end

    case 'c'

        eta = varargin{1};
        [Sxx, Sxy, Syy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 's');
        [Dxx, Dxy, Dyy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'd');
        Kxx = Dxx + eta*Sxx;
        Kyy = Dyy + eta*Syy;
        Kxy = Dxy + eta*Sxy;

        if ( nargout == 1 )
            % Interleave
            K = zeros(2*Nt, 2*Ns);
            K(1:2:end, 1:2:end) = Kxx;
            K(1:2:end, 2:2:end) = Kxy;
            K(2:2:end, 1:2:end) = Kxy;
            K(2:2:end, 2:2:end) = Kyy;
            varargout = {K};
        else
            varargout = {Kxx, Kxy, Kyy};
        end
    case 'cpres'

        eta = varargin{1};
        [Sx, Sy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'spres');
        [Dx, Dy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'dpres');
        Kx = Dx + eta*Sx;
        Ky = Dy + eta*Sy;    

        if ( nargout == 1 )
            % Interleave
            K = zeros(Nt, 2*Ns);
            K(:,1:2:end) = Kx;
            K(:,2:2:end) = Ky;
            varargout = {K};
        else
            varargout = {Kx, Ky};
        end

    case 'ctrac'

        eta = varargin{1};
        [Sxx, Sxy, Syy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'strac');
        [Dxx, Dxy, Dyy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'dtrac');
        Kxx = Dxx + eta*Sxx;
        Kxy = Dxy + eta*Sxy;
        Kyy = Dyy + eta*Syy;
        
        if ( nargout == 1 )
            % Interleave
            K = zeros(2*Nt, 2*Ns);
            K(1:2:end, 1:2:end) = Kxx;
            K(1:2:end, 2:2:end) = Kxy;
            K(2:2:end, 1:2:end) = Kxy;
            K(2:2:end, 2:2:end) = Kyy;
            varargout = {K};
        else
            varargout = {Kxx, Kxy, Kyy};
        end

        

    otherwise
        error('Unknown Stokes kernel type ''%s''.', type);

end

end
