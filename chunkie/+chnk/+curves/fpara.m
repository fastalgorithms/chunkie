function [r,d,d2] = fpara(t,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
%   a*(t-b)^2
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ts = size(t);
    tt = t(:);
    
    rx = tt;
    ry = a*(tt-b).^2;
    
    dx = 1*ones(size(rx));
    dy = 2*a*(tt-b);
    
    d2x= zeros(size(rx));
    d2y= 2*a(ones(size(d2x)));
    
    r = [rx.';ry.'];
    d = [dx.';dy.'];
    d2= [d2x.';d2y.'];
    
    r = reshape(r,[2,ts]);
    d = reshape(d,[2,ts]);
    d2= reshape(d2,[2,ts]);
    
end

