function obj = reverse(obj)
%REVERSE reverses the orientation of the trapper object
%

obj.r = obj.r(:,end:-1:1);
obj.d = -obj.d(:,end:-1:1);
obj.d2 = obj.d2(:,end:-1:1);

end