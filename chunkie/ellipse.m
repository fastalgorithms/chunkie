function [r, d, d2] = ellipse(t)

x = 3.*cos(t);
y = sin(t);

x1 = -3.*sin(t);
y1 = cos(t);

x2 = -3.*cos(t);
y2 = -sin(t);

r = [(x(:)).'; (y(:)).'];
d = [(x1(:)).'; (y1(:)).'];
d2 = [(x2(:)).'; (y2(:)).'];

end