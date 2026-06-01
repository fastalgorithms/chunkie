function cg = reflect(cg,theta,r0,r1)
%REFLECT reflect chunkgraph object across line specified by angle
%
% Syntax 
%   cg = cg.reflect(theta,r0,r1);
%
% will shift the chunkgraph geometry by -r0, reflect it across the line 
% spanned by [cos(theta); sin(theta)] and shift back by +r1, i.e. the 
% coordinates will be changed by 
%
% r <- M*(r-r0) + r1
%
% where M is the reflection matrix [cos(2*theta) sin(2*theta); sin(2*theta)
% -cos(2*theta)]
%
% Input:
%   theta - float, angle made by reflection axis in radians. if empty, 
%              theta=0
%   r0 - length 2 array, center of reflection if empty, r0=[0;0]
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

assert(isreal(theta),"reflect only supports real angles");

for j = 1:length(cg.echnks)
    cg.echnks(j) = cg.echnks(j).reflect(theta,r0,r1);
end

c = cos(2*theta);
s = sin(2*theta);
refmat = [c s; s -c];

cg.verts = refmat*(cg.verts-r0) + r1;

% region stuff needs recomp
cg.vstruc = procverts(cg);
cg.regions = findregions(cg);
