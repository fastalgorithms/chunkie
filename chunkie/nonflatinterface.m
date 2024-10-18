function [r,d,d2] = nonflatinterface(t,a,b,c,d)
   xs = t;
   xp = ones(size(t));
   xpp = zeros(size(t));
   
   
    ys = d*exp(-a*t.^2/2).*sin(b*t+c);
    yp = d*exp(-a*t.^2/2).*(b*cos(b*t+c) - a*t.*sin(b*t+c));
    ypp = d*exp(-a*t.^2/2).*(-2*a*b*t.*cos(b*t+c) + (a*a*t.^2-a-b^2).*sin(b*t+c));
    
    r = [(xs(:)).'; (ys(:)).'];
    d = [(xp(:)).'; (yp(:)).'];
    d2 = [(xpp(:)).'; (ypp(:)).'];
end