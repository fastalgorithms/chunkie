%GUIDE_CHUNKERS
%
% This script complements the chunkie guide section on chunker objects.
% It shows the most common methods for building and working with chunkers.
%

addpaths_loc();
rng(1234);

%%%%%%%%%%%%% circle
% START CIRCLE
% chunk up circle

rad = 2; ctr = [1.0;-0.5];
circfun = @(t) ctr + rad*[cos(t(:).');sin(t(:).')];

chnkr1 = chunkerfunc(circfun);

% plot curve, nodes, and normals

figure(1)
clf
plot(chnkr1,'b-x')
hold on
quiver(chnkr1,'r')
axis equal tight
% END CIRCLE

saveas(figure(1),"guide_chunkers_circle.png");

%%%%%%%%%%%%% more examples 
% START MORE PARAMS
% curvebymode lets you specify a star shaped domain by 
% a cosine/sine series for the radius

modes = randn(11,1); modes(1) = 1.1*sum(abs(modes(2:end))); 
ctr = [1.0;-0.5];
chnkr2 = chunkerfunc(@(t) chnk.curves.bymode(t,modes,ctr));

figure(2)
clf
plot(chnkr2,'b-x')
hold on
quiver(chnkr2,'r')
axis equal tight

% classic starfish domain

narms = 5; 
amp = 0.5;
chnkr3 = chunkerfunc(@(t) starfish(t,narms,amp));

figure(3)
clf
plot(chnkr3,'b-x')
hold on
quiver(chnkr3,'r')
axis equal tight
% END MORE PARAMS

saveas(figure(2),"guide_chunkers_bymode.png");
saveas(figure(3),"guide_chunkers_starfish.png");

%%%%%%%%%%%%% rounded polygons
% START ROUNDED POLY
% chunkerpoly provides a chunk discretization of a rounded polygon

verts = chnk.demo.barbell(2.0,2.0,1.0,1.0); % vertices of a barbell
chnkr4 = chunkerpoly(verts);

figure(4)
clf
plot(chnkr4,'b-x')
hold on
quiver(chnkr4,'r')
axis equal tight
% END ROUNDED POLY

saveas(figure(4),"guide_chunkers_barbell.png");

%%%%%%%%%%%%%% working with chunkers

% START SHIFT AND REVERSE
% make a copy of the random mode domain
chnkr5 = chnkr2;

% rotate it using move 
theta = pi/4; 
chnkr5 = chnkr5.move([],[],theta);

% make a copy of the circle domain and shift it 
chnkr6 = chnkr1;
r1 = [-0.3;0.4];
chnkr6 = chnkr6.move([],r1);
% reverse the orientation
chnkr6 = chnkr6.reverse();

% merge these two curves into one domain
chnkr7 = merge([chnkr5,chnkr6]);

figure(5); clf
plot(chnkr7,'b-x'); hold on; quiver(chnkr7,'r')
axis equal tight
% END SHIFT AND REVERSE

saveas(figure(5),"guide_chunkers_shiftandreverse.png");

% START INTERIOR

% create a grid of points 
mins = min(chnkr7); maxs = max(chnkr7);
x1 = linspace(mins(1),maxs(1)); y1 = linspace(mins(2),maxs(2));
[xx,yy] = meshgrid(x1,y1);
pts = [xx(:).'; yy(:).'];

% find the points inside the merged domain
in = chunkerinterior(chnkr7,pts);

zz = nan(size(xx));
zz(in) = 1;

figure(6)
h = pcolor(xx,yy,zz); set(h,'EdgeColor','none');
axis equal tight
% END INTERIOR

saveas(figure(6),"guide_chunkers_interior.png");

