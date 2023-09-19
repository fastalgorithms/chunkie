function obj = move(obj,r0,r1,trotat,scale)

    
rsize   = size(obj.r);

rotmat = [cos(trotat),-sin(trotat);sin(trotat),cos(trotat)];

rnew = scale*rotmat*(obj.r(:,:)-r0) + r1;
obj.r = reshape(rnew,rsize);

dnew = scale*rotmat*(obj.d(:,:));
obj.d = reshape(dnew,rsize);

d2new = scale*rotmat*(obj.d2(:,:));
obj.d2 = reshape(d2new,rsize);

nnew = rotmat*(obj.n(:,:));
obj.n = reshape(nnew,rsize);

obj.wts = weights(obj);
 
end