function [r,d,d2] = fsine(t,a,b,c)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
%   a*sin(b*t+c)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ts = size(t);
    tt = t(:);
    
    rx = tt;
    ry = a*sin(b*tt+c);
    
    dx = 1*ones(size(rx));
    dy = a*b*cos(b*tt+c);
    
    d2x= zeros(size(rx));
    d2y=-a*b*b*sin(b*tt+c);
    
    r = [rx.';ry.'];
    d = [dx.';dy.'];
    d2= [d2x.';d2y.'];
    
    r = reshape(r,[2,ts]);
    d = reshape(d,[2,ts]);
    d2= reshape(d2,[2,ts]);
    
end

