 




% define the geometry of the circle 
function [r, d, d2]= circle(t) 

x = cos(t);
y = sin(t);

dx = -sin(t);
dy = cos(t);

dxx = cos(t);
dyy = sin(t);

r = [(x(:)).';(y(:)).'];
d = [(dx(:)).'; (dy(:)).'];
d2 = [(dxx(:)).'; (dyy(:)).'];


end