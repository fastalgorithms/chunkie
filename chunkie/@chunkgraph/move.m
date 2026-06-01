function obj = move(obj, r0, r1, trotat, scale)
%MOVE update the position of chunkers in the chunkgraph by 
% translation, rotation, and scaling.
%
% Syntax 
%   chnkgraph = chnkgraph.move(r0,r1,trotat,scale)
%
% will shift the chunker geometry by -r0, rotate it by angle trotat, 
% scale the points by the factor scale, and shift back by +r1, i.e. the
% coordinates will be changed by 
%
% r <- s*M*(r-r0) + r1
%
% where M is the rotation matrix for trotat radians and s is the scale
%
% Input:
%   r0 - length 2 array, center of rotation. if empty, r0=[0;0]
%   r1 - length 2 array, new center. if empty, r1=[0;0]
%   trotat - float, angle of rotation in radians. if empty, trotat=0
%   scale - float, scaling factor. ifempty, scale = 1
% 
% Output:
%   chnkr - chunker object with positions, derivatives, normals etc updated
%   according to the formula above. 
%

if nargin < 2 || isempty(r0)
    r0 = [0;0];
end
if nargin < 3 || isempty(r1)
    r1 = [0;0];
end
if nargin < 4 || isempty(trotat)
    trotat = 0;
end
if nargin < 5 || isempty(scale)
    scale = 1;
end

for i = 1:length(obj.echnks)
    obj.echnks(i) = move(obj.echnks(i), r0,r1,trotat,scale);
end

vsize   = size(obj.verts);

rotmat = [cos(trotat),-sin(trotat);sin(trotat),cos(trotat)];

vertsnew = scale*rotmat*(obj.verts(:,:)-r0) + r1;
obj.verts = reshape(vertsnew, vsize);
 


end
