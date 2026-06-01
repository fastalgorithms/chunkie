function cg = rotate(cg,theta,r0,r1)
%ROTATE rotate chunkgraph object by specified angle
%
% Syntax 
%   cg = cg.rotate(theta,r0,r1);
%
% will shift the chunkgraph geometry by -r0, rotate it by angle trotat, 
% and shift back by +r1, i.e. the coordinates will be changed by 
%
% r <- M*(r-r0) + r1
%
% where M is the rotation matrix for theta radians
%
% Input:
%   theta - float, angle of rotation in radians. if empty, trotat=0
%   r0 - length 2 array, center of rotation. if empty, r0=[0;0]
%   r1 - length 2 array, new center. if empty, r1=[0;0]
% 
% Output:
%   cg - chunkgraph object with positions, derivatives, normals etc updated
%   according to the formula above. 
%

if nargin < 3 || isempty(r0)
    r0 = [0;0];
end
if nargin < 4 || isempty(r1)
    r1 = [0;0];
end

assert(isreal(theta),"rotate only supports real angles");

for j = 1:length(cg.echnks)
    cg.echnks(j) = cg.echnks(j).rotate(theta,r0,r1);
end

c = cos(theta);
s = sin(theta);
rotmat = [c -s; s c];

cg.verts = rotmat*(cg.verts-r0) + r1;