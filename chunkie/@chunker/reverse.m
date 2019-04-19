function obj = reverse(obj)
%REVERSE reverses the orientation of the chunker object
%

k = obj.k;
obj.r = obj.r(:,k:-1:1,:);
obj.d = -obj.d(:,k:-1:1,:);
obj.d2 = obj.d2(:,k:-1:1,:);
obj.adj = obj.adj([2 1],:);

end