function varargout = kern(mu, srcinfo, targinfo, type, varargin)
%CHNK.STOK2D.KERN   Evaluate the Stokes kernel.

rx = targinfo.r(1,:).' - srcinfo.r(1,:);
ry = targinfo.r(2,:).' - srcinfo.r(2,:);
r2 = rx.^2 + ry.^2;
[Nt, Ns] = size(rx);

switch lower(type)

    case {'svel', 's', 'single'}

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
    
    case {'sgrad', 'sg'}
        r4 = r2.^2;
        r4inv = 1./r4./(4*pi*mu);
        Ku1x1 = rx.*(ry.^2 - rx.^2).*r4inv;
        Ku1x2 = ry.*(ry.^2 - rx.^2).*r4inv;

        Ku1y1 = (-3.*rx.^2.*ry - ry.^3).*r4inv;
        Ku1y2 = rx.*(rx.^2 - ry.^2).*r4inv;

        Ku2x1 = ry.*(ry.^2 - rx.^2).*r4inv;
        Ku2x2 = (-3.*rx.*ry.^2 - rx.^3).*r4inv;

        Ku2y1 = rx.*(rx.^2 - ry.^2).*r4inv;
        Ku2y2 = ry.*(rx.^2 - ry.^2).*r4inv;

        if ( nargout == 1 )
            % Interleave
            K = zeros(4*Nt, 2*Ns);
            K(1:4:end, 1:2:end) = Ku1x1;
            K(1:4:end, 2:2:end) = Ku1x2;
            K(2:4:end, 1:2:end) = Ku1y1;
            K(2:4:end, 2:2:end) = Ku1y2;

            K(3:4:end, 1:2:end) = Ku2x1;
            K(3:4:end, 2:2:end) = Ku2x2;
            K(4:4:end, 1:2:end) = Ku2y1;
            K(4:4:end, 2:2:end) = Ku2y2; 

            varargout = {K};
        else
            varargout = {Ku1x1, Ku1x2, Ku1y1, Ku1y2, Ku2x1, Ku2x2, Ku2y1, Ku2y2};
        end



    case {'dvel', 'd', 'double'}

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
            nx_t.*nx_s ./ r2);
        
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

    case {'dgrad', 'dg'}
        r4 = r2.^2;
        r6 = r2.^3;
        nx = srcinfo.n(1,:);
        ny = srcinfo.n(2,:);
        rn = rx.*nx + ry.*ny;
        overpi = 1/pi;

        r4inv = 1./r4;
        r6inv = 1./r6;
        Ku1x1 = (rn.*rx.*r4inv + rx.*rn.*r4inv + rx.^2.*nx.*r4inv - 4*rx.^3.*rn.*r6inv)*overpi;
        Ku1x2 = (rn.*ry.*r4inv + rx.*ry.*nx.*r4inv - 4*rx.^2.*ry.*rn.*r6inv)*overpi;

        Ku1y1 = (rx.^2.*ny.*r4inv - 4*rx.^2.*ry.*rn.*r6inv)*overpi;
        Ku1y2 = (rx.*rn.*r4inv + rx.*ry.*ny.*r4inv - 4*rx.*ry.^2.*rn.*r6inv)*overpi;

        Ku2x1 = (ry.*rn.*r4inv + rx.*ry.*nx.*r4inv - 4*rx.^2.*ry.*rn.*r6inv)*overpi;
        Ku2x2 = (ry.^2.*nx.*r4inv - 4*rx.*ry.^2.*rn.*r6inv)*overpi;

        Ku2y1 = (rn.*rx.*r4inv + rx.*ry.*ny.*r4inv - 4*rx.*ry.^2.*rn.*r6inv)*overpi;
        Ku2y2 = (rn.*ry.*r4inv + ry.*rn.*r4inv + ry.^2.*ny.*r4inv - 4*ry.^3.*rn.*r6inv)*overpi;

        if ( nargout == 1 )
            % Interleave
            K = zeros(4*Nt, 2*Ns);
            K(1:4:end, 1:2:end) = Ku1x1;
            K(1:4:end, 2:2:end) = Ku1x2;
            K(2:4:end, 1:2:end) = Ku1y1;
            K(2:4:end, 2:2:end) = Ku1y2;

            K(3:4:end, 1:2:end) = Ku2x1;
            K(3:4:end, 2:2:end) = Ku2x2;
            K(4:4:end, 1:2:end) = Ku2y1;
            K(4:4:end, 2:2:end) = Ku2y2; 


            varargout = {K};
        else
            varargout = {Ku1x1, Ku1x2, Ku1y1, Ku1y2, Ku2x1, Ku2x2, Ku2y1, Ku2y2};
        end


    case {'c', 'combined'}

        coefs = varargin{1};
        [Sxx, Sxy, Syy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 's');
        [Dxx, Dxy, Dyy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'd');
        Kxx = coefs(1)*Dxx + coefs(2)*Sxx;
        Kyy = coefs(1)*Dyy + coefs(2)*Syy;
        Kxy = coefs(1)*Dxy + coefs(2)*Sxy;

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

        coefs = varargin{1};
        [Sx, Sy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'spres');
        [Dx, Dy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'dpres');
        Kx = coefs(1)*Dx + coefs(2)*Sx;
        Ky = coefs(1)*Dy + coefs(2)*Sy;    

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

        coefs = varargin{1};
        [Sxx, Sxy, Syy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'strac');
        [Dxx, Dxy, Dyy] = chnk.stok2d.kern(mu, srcinfo, targinfo, 'dtrac');
        Kxx = coefs(1)*Dxx + coefs(2)*Sxx;
        Kxy = coefs(1)*Dxy + coefs(2)*Sxy;
        Kyy = coefs(1)*Dyy + coefs(2)*Syy;
        
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

    case {'cgrad', 'cg'}
        coefs = varargin{1};
        [Su1x1, Su1x2, Su1y1, Su1y2, Su2x1, Su2x2, Su2y1, Su2y2] = ...
           chnk.stok2d.kern(mu, srcinfo, targinfo, 'sgrad');
        [Du1x1, Du1x2, Du1y1, Du1y2, Du2x1, Du2x2, Du2y1, Du2y2] = ...
           chnk.stok2d.kern(mu, srcinfo, targinfo, 'sgrad');

        Ku1x1 = coefs(1)*Du1x1 + coefs(2)*Su1x1;
        Ku1x2 = coefs(1)*Du1x2 + coefs(2)*Su1x2;

        Ku1y1 = coefs(1)*Du1y1 + coefs(2)*Su1y1;
        Ku1y2 = coefs(1)*Du1y2 + coefs(2)*Su1y2;

        Ku2x1 = coefs(1)*Du2x1 + coefs(2)*Su2x1;
        Ku2x2 = coefs(1)*Du2x2 + coefs(2)*Su2x2;

        Ku2y1 = coefs(1)*Du2y1 + coefs(2)*Su2y1;
        Ku2y2 = coefs(1)*Du2y2 + coefs(2)*Su2y2;

        if ( nargout == 1 )
            % Interleave
            K = zeros(4*Nt, 2*Ns);
            K(1:4:end, 1:2:end) = Ku1x1;
            K(1:4:end, 2:2:end) = Ku1x2;
            K(2:4:end, 1:2:end) = Ku1y1;
            K(2:4:end, 2:2:end) = Ku1y2;

            K(3:4:end, 1:2:end) = Ku2x1;
            K(3:4:end, 2:2:end) = Ku2x2;
            K(4:4:end, 1:2:end) = Ku2y1;
            K(4:4:end, 2:2:end) = Ku2y2; 

            varargout = {K};
        else
            varargout = {Ku1x1, Ku1x2, Ku1y1, Ku1y2, Ku2x1, Ku2x2, Ku2y1, Ku2y2};
        end        

    otherwise
        error('Unknown Stokes kernel type ''%s''.', type);

end

end
