function verts = barbell(s1,s2,hb,lb)
%BARBELL
%
% Input (all optional)
%   s1 - side length of left square (4.0)
%   s2 - side length of right square (3.0)
%   hb - height of bridge between squares (1.0)
%   lb - length of bridge (2.5)
%

if nargin < 1; s1 = 2.0; end
if nargin < 2; s2 = 2.0; end
if nargin < 3; hb = 0.5; end
if nargin < 4; lb = 1.0; end

verts_sqr1 = s1/2.0*[1, -1, -1, 1; 1, 1, -1, -1];
verts_sqr2 = s2/2.0*[-1, 1, 1, -1; -1, -1, 1, 1];

wbrdg = hb/2.0;

ctr_sqr1 = [-(s1/2.0+lb/2.0); 0.0];
ctr_sqr2 = [(s2/2.0+lb/2.0); 0.0];

verts_sqr1 = verts_sqr1 + repmat(ctr_sqr1,1,size(verts_sqr1,2));
verts_sqr2 = verts_sqr2 + repmat(ctr_sqr2,1,size(verts_sqr2,2));

verts1 = [verts_sqr1(1,end), verts_sqr2(1,1); -wbrdg, -wbrdg];
verts2 = [verts_sqr2(1,end), verts_sqr1(1,1); wbrdg, wbrdg];

verts = [verts_sqr1, verts1, verts_sqr2, verts2];